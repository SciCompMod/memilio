/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Henrik Zunker
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef MIO_ODE_SEIR_WORLD_APPLICATION_H
#define MIO_ODE_SEIR_WORLD_APPLICATION_H

#include "ode_seir_runtime_scenario.h"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mio::world_application
{

// Worldwide application. Populations, initial states and the mobility matrix are read
// from the binary file written by hpc_paper/world/preprocess_world.py, the epidemiological
// parameters are those of the runtime scenario. Enabled by WORLD_INPUT=<file>.
inline const char* input_path()
{
    return std::getenv("WORLD_INPUT");
}

inline bool enabled()
{
    const char* path = input_path();
    return path != nullptr && *path != '\0';
}

inline runtime_scenario::Schedule schedule()
{
    return {runtime_scenario::positive_environment("WORLD_DAYS", 365, 3660),
            runtime_scenario::positive_environment("RUNTIME_HALF_DAY_STEPS", 32, 1024)};
}

inline std::string output_directory()
{
    const char* directory = std::getenv("WORLD_OUTPUT_DIR");
    return directory != nullptr && *directory != '\0' ? directory : ".";
}

inline runtime_scenario::Inputs load_inputs()
{
    std::ifstream stream(input_path(), std::ios::binary);
    std::int32_t shape[2] = {0, 0};
    stream.read(reinterpret_cast<char*>(shape), sizeof(shape));
    if (!stream) {
        throw std::runtime_error("Cannot read the world input header.");
    }
    // The constructor checks the shape and sets the epidemiological parameters.
    runtime_scenario::Inputs inputs(shape[0], shape[1]);
    for (auto* values : {&inputs.population, &inputs.initial, &inputs.h}) {
        stream.read(reinterpret_cast<char*>(values->data()),
                    static_cast<std::streamsize>(values->size() * sizeof(double)));
    }
    if (!stream || stream.peek() != std::ifstream::traits_type::eof()) {
        throw std::runtime_error("The world input size does not match its header.");
    }
    const size_t patches = static_cast<size_t>(inputs.p);
    for (size_t origin = 0; origin < patches; ++origin) {
        for (size_t destination = 0; destination < patches; ++destination) {
            inputs.ht[destination * patches + origin] = inputs.h[origin * patches + destination];
        }
    }
    return inputs;
}

// Infectious residents per patch, summed over the age groups, at the end of each day.
class DailyInfectious
{
public:
    explicit DailyInfectious(const runtime_scenario::Inputs& inputs)
        : m_inputs(inputs)
    {
    }

    void record(const std::vector<double>& residents)
    {
        for (int patch = 0; patch < m_inputs.p; ++patch) {
            double infectious = 0.0;
            for (int group = 0; group < m_inputs.g; ++group) {
                infectious += residents[m_inputs.index(patch, group, 2)];
            }
            m_values.push_back(static_cast<float>(infectious));
        }
        ++m_records;
    }

    // Binary file: int32 records, int32 patches, float32 values[records * patches].
    void write(const std::string& model, runtime_scenario::Schedule time, double seconds) const
    {
        const std::string path = output_directory() + "/" + model + "_daily_infectious.bin";
        std::ofstream stream(path, std::ios::binary);
        const std::int32_t shape[2] = {m_records, m_inputs.p};
        stream.write(reinterpret_cast<const char*>(shape), sizeof(shape));
        stream.write(reinterpret_cast<const char*>(m_values.data()),
                     static_cast<std::streamsize>(m_values.size() * sizeof(float)));
        if (!stream) {
            throw std::runtime_error("Cannot write " + path + ".");
        }
        std::cout << "World " << model << ": patches=" << m_inputs.p << " age_groups=" << m_inputs.g
                  << " days=" << time.days << " dt=" << time.dt() << " seconds=" << seconds << " output=" << path
                  << '\n';
    }

private:
    const runtime_scenario::Inputs& m_inputs;
    std::vector<float> m_values;
    std::int32_t m_records = 0;
};

// Simulates the given number of days and records the residents after each day. The
// measured time contains only the integration, not the daily download of the states.
template <class AdvanceDay, class Residents>
void simulate(const runtime_scenario::Inputs& inputs, const std::string& model, runtime_scenario::Schedule time,
              AdvanceDay advance_day, Residents residents)
{
    DailyInfectious output(inputs);
    output.record(inputs.initial);
    std::vector<double> state = inputs.initial;
    double seconds            = 0.0;
    for (int day = 0; day < time.days; ++day) {
        const auto start = std::chrono::steady_clock::now();
        advance_day();
        seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        state = residents();
        output.record(state);
    }
    runtime_scenario::check_population(inputs, state);
    output.write(model, time, seconds);
}

} // namespace mio::world_application

#endif // MIO_ODE_SEIR_WORLD_APPLICATION_H
