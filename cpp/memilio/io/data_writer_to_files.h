/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: René Schmieding
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
#ifndef MIO_IO_DATA_WRITER_TO_FILES_H
#define MIO_IO_DATA_WRITER_TO_FILES_H

#include "memilio/config.h" // IWYU pragma: keep
#include <cstddef>

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/json_serializer.h"
#include "memilio/utils/logging.h"

namespace mio
{

/**
 * @brief Writer that writes records from loggers directly to files.
 * Uses exactly one file per Logger. 
 * Prefers the stream operator "<<" to write a record if available, falls back to json serialization. 
 */
template <class... Loggers>
struct DataWriterToFiles {
    /**
     * @brief Struct to hold the output streams and a json writer.
     * Provides a constructor for generic ostreams, and for file streams (taking file paths as strings).
     * The json writer uses minimal white space to hopefully reduce the cost of writing (especially for files).
     */
    struct Data {
        /// @brief Number of Loggers
        static constexpr size_t N = sizeof...(Loggers);

        /// @brief Take over an existing set of std::ofstream%s.
        Data(std::array<std::ofstream, N> log_streams)
            : outputs(std::move(log_streams))
        {
            Json::StreamWriterBuilder builder;
            builder["indentation"] = "";
            writer.reset(builder.newStreamWriter());
        }

        /// @brief Open new std::ofstream%s from the given file paths.
        Data(std::array<std::string, N> filepaths)
            : Data([](auto&& paths) {
                std::array<std::ofstream, N> files;
                for (size_t i = 0; i < N; i++) {
                    files[i].open(paths[i]);
                    if (!files[i].good()) {
                        log_critical("Error opening file \"{}\" in DataWriterToOutputStreams.", paths[i]);
                    }
                }
                return files;
            }(filepaths))
        {
        }

        /// @brief Test whether all output streams are good.
        bool good() const
        {
            return std::ranges::all_of(outputs, [](auto&& o) {
                return o.good();
            });
        }

        std::array<std::ofstream, N> outputs;
        std::unique_ptr<Json::StreamWriter> writer;
    };

    /**
     * @brief Adds a new record for a given log result t to data.
     * The parameter Logger is used to determine the type of the record t, as well as the data index at which
     * the record should be added to.
     * @param[in] t The result of Logger::log.
     * @param[in,out] data An instance of Data to add the record to.
     * @tparam Logger The type of the logger used to record t.
     */
    template <class Logger>
    static void add_record(const typename Logger::Type& t, Data& data)
    {
        auto& out = data.outputs[mio::index_of_type_v<Logger, Loggers...>];
        if constexpr (requires { std::cout << t; }) {
            out << t;
        }
        else {
            data.writer->write(serialize_json(t).value(), &out);
        }
        out << "\n";
    }
};

} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // MIO_IO_DATA_WRITER_TO_FILES_H
