/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Rene Schmieding
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
#ifndef BENCH_CORE_H_
#define BENCH_CORE_H_

#include "bench_model.h"

#include "memilio/utils/logging.h"
#include "memilio/math/adapt_rk.h"
#include "memilio/math/stepper_wrapper.h"

#include <benchmark/benchmark.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>

namespace mio
{
namespace benchmark
{
    namespace detail
    {
        inline void reset_time(double& t, double& dt, const double& t_init, const double& dt_init) {
            t  = t_init;
            dt = dt_init;
        }
        /**
         * @brief Basic file handler to read parameters from a config file.
         */
        class ConfigReader
        {
        public:
            ConfigReader(char const *filepath) : m_file(filepath), m_filename(filepath) {
                if (!m_file.good()) {
                    mio::log(mio::LogLevel::critical,
                             " Could not open config file \'{:s}\'.\n",
                             m_filename);
                    abort();
                }
            }
            ~ConfigReader() {
                m_file.close();
            }
            /**
             * @brief Converts a line from the provided config file into data
             */
            template <class dtype>
            void read(dtype& data) {
                std::string reader;
                getline<dtype>(reader);
                // use stringstream to convert to various data types
                std::stringstream parser;
                parser << reader;
                parser >> data;
            }
            /**
             * @brief Converts a line from the provided config file into data
             */
            void read(double& data) {
                std::string reader;
                getline<double>(reader);
                if (reader.compare("max") == 0) {
                    data = std::numeric_limits<double>::max();
                }
                else if (reader.compare("min") == 0) {
                    data = std::numeric_limits<double>::min();
                }
                else {
                    std::stringstream parser;
                    parser << reader;
                    parser >> data;
                }
            }
            /**
             * @brief Converts a line from the provided config file into data
             */
            void read(std::vector<double>& data) {
                std::string reader;
                std::stringstream parser;
                getline<double>(reader);
                parser << reader;
                for (double& d : data) {
                    parser >> d;
                }
            }
        private:
            std::ifstream m_file;
            char const *m_filename;
            /**
             * @brief calls std::getline with internal filestream and logs / aborts on failure
             */
            template <class dtype>
            void getline(std::string& reader) {
                if (!std::getline(m_file, reader)) {
                    mio::log(mio::LogLevel::critical,
                             " Unexpected end of file in \'{:s}\'.\n Expected a value of type \'{:s}\'.\n",
                             m_filename, typeid(dtype).name());
                    abort();
                }
            }
        }; // class ConfigReader
    } // namespace detail

    template <class Model=mio::benchmark::model::SecirAgeres, size_t num_agegroups=1>
    struct default_simulation_init
    {
        typename Model::type operator () (
            double& abs_tol, double& rel_tol, double& dt_min, double& dt_max,
            double& t, double& t_max, double& dt
        ) {
            t       = 0;
            t_max   = 100;
            dt      = 0.5;
            abs_tol = 1e-10;
            rel_tol = 1e-5;
            dt_min  = std::numeric_limits<double>::min();
            dt_max  = std::numeric_limits<double>::max();
            return Model()(num_agegroups);
        }
    };

    template <char const *filepath, class Model=mio::benchmark::model::SecirAgeres>
    struct simulation_file_init
    {
        typename Model::type operator () (
            double& abs_tol, double& rel_tol, double& dt_min, double& dt_max,
            double& t, double& t_max, double& dt
        ) {
            // read configuration file
            mio::benchmark::detail::ConfigReader config(filepath);
            size_t num_agegroups;
            config.read(num_agegroups);
            config.read(t);
            config.read(t_max);
            config.read(dt);
            config.read(abs_tol);
            config.read(rel_tol);
            config.read(dt_min);
            config.read(dt_max);
            return Model()(num_agegroups);
        }
    };

    template <class Integrator, class Initializer=default_simulation_init<>>
    void simulation(::benchmark::State& state)
    {
        // suppress non-critical messages
        mio::set_log_level(mio::LogLevel::critical);
        // benchmark setup
        double abs_tol, rel_tol, dt_min, dt_max, t, t_max, dt;

        auto model = Initializer()(abs_tol, rel_tol, dt_min, dt_max, t, t_max, dt);

        for (auto _ : state) {
            // This code gets timed
            std::shared_ptr<mio::IntegratorCore> I = std::make_shared<Integrator>(abs_tol, rel_tol, dt_min, dt_max);
            simulate(t, t_max, dt, model, I);
        }
    }

    template <class Model=mio::benchmark::model::SecirAgeres>
    struct default_integrator_step_init
    {
        typename Model::type operator () (
            double& abs_tol, double& rel_tol, double& dt_min, double& dt_max,
            double& t_init, double& dt_init,
            mio::DerivFunction& f, Eigen::VectorXd& yt, Eigen::VectorXd& ytp1
        ) {
            typename Model::type model = Model()(1);
            t_init  = 50;
            dt_init = 1;
            abs_tol = 1e-10;
            rel_tol = 1e-5;
            dt_min  = std::numeric_limits<double>::min();
            dt_max  = std::numeric_limits<double>::max();
            f       = [model](Eigen::Ref<const Eigen::VectorXd> y_f, double t_f, Eigen::Ref<Eigen::VectorXd> dydt_f) { model.eval_right_hand_side(y_f, y_f, t_f, dydt_f); };
            yt      = Eigen::VectorXd::Zero(model.get_initial_values().size());
            ytp1    = Eigen::VectorXd::Zero(yt.size());
            // set values for yt
            std::vector<double> yt_vals({6377.873644, 35.249156, 30.029611, 182.145865, 66.153059, 79.530621, 3069.383604, 159.634440});
            assert(yt.size() <= yt_vals.size());
            for (int i = 0; i < yt.size(); i++) {
                yt[i] = yt_vals[i];
            }
            return model;
        }
    };

    template <char const *filepath, class Model=mio::benchmark::model::SecirAgeres>
    struct integrator_step_file_init
    {
        typename Model::type operator () (
            double& abs_tol, double& rel_tol, double& dt_min, double& dt_max,
            double& t_init, double& dt_init,
            mio::DerivFunction& f, Eigen::VectorXd& yt, Eigen::VectorXd& ytp1
        ) {
            // read configuration file
            mio::benchmark::detail::ConfigReader config(filepath);
            size_t num_agegroups;
            config.read(num_agegroups);
            config.read(t_init);
            config.read(dt_init);
            config.read(abs_tol);
            config.read(rel_tol);
            config.read(dt_min);
            config.read(dt_max);
            // create model and read initial yt values
            typename Model::type model = Model()(num_agegroups);
            std::vector<double> yt_vals(model.get_initial_values().size());
            config.read(yt_vals);
            // initialize remaining components
            f     = [model](Eigen::Ref<const Eigen::VectorXd> y_f, double t_f, Eigen::Ref<Eigen::VectorXd> dydt_f) { model.eval_right_hand_side(y_f, y_f, t_f, dydt_f); };
            yt    = Eigen::VectorXd::Zero(yt_vals.size());
            ytp1  = Eigen::VectorXd::Zero(yt_vals.size());
            for (int i = 0; i < yt.size(); i++) {
                yt[i] = yt_vals[i];
            }
            return model;
        }
    };

    template <class Integrator, class Initializer=default_integrator_step_init<>>
    void integrator_step(::benchmark::State& state)
    {
        // suppress non-critical messages
        mio::set_log_level(mio::LogLevel::critical);
        // benchmark setup
        double abs_tol, rel_tol, dt_min, dt_max, t_init, dt_init;
        mio::DerivFunction f;
        Eigen::VectorXd yt, ytp1;

        auto model = Initializer()(abs_tol, rel_tol, dt_min, dt_max, t_init, dt_init, f, yt, ytp1);

        auto I = Integrator(abs_tol, rel_tol, dt_min, dt_max);

        double t, dt;
        for (auto _ : state) {
            // This code gets timed
            mio::benchmark::detail::reset_time(t, dt, t_init, dt_init);
            I.step(f, yt, t, dt, ytp1);
        }
    }
} // namespace benchmark

} // namespace mio

#endif