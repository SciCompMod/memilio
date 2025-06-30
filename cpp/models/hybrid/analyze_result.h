/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker
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
#ifndef MIO_HYBRID_ANALYZE_RESULT_H
#define MIO_HYBRID_ANALYZE_RESULT_H

#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"

namespace mio
{

namespace hybrid
{

/**
 * @brief This function merges two TimeSeries by copying their time points and values to a new TimeSeries in the correct order.
 * If both TimeSeries have values for the same time point, their values are either added or only one value is taken.
 * @param[in] ts1 First TimeSeries.
 * @param[in] ts2 Second TimeSeries.
 * @param[in] add_values Boolean specifying whether the values should be added if both TimeSeries contain the same time point. If false, the value of just the first TimeSeries is taken.
 * @tparam FP A floating point type.
 * @return A TimeSeries containing all time points and values from both input TimeSeries.
 */
template <class FP>
TimeSeries<FP> merge_time_series(TimeSeries<FP>& ts1, TimeSeries<FP>& ts2, bool add_values = false)
{
    TimeSeries<FP> merged_ts(ts1.get_num_elements());
    if (ts1.get_num_elements() != ts2.get_num_elements()) {
        log_error("TimeSeries have a different number of elements.");
    }
    else {
        Eigen::Index t1_iterator = 0;
        Eigen::Index t2_iterator = 0;
        bool t1_finished         = false;
        bool t2_finished         = false;
        while (!t1_finished || !t2_finished) {
            if (!t1_finished) {
                if (ts1.get_time(t1_iterator) < ts2.get_time(t2_iterator) ||
                    t2_finished) { // Current time point of first TimeSeries is smaller than current time point of second TimeSeries or second TimeSeries has already been copied entirely
                    merged_ts.add_time_point(ts1.get_time(t1_iterator), ts1.get_value(t1_iterator));
                    t1_iterator += 1;
                }
                else if (!t2_finished && ts1.get_time(t1_iterator) ==
                                             ts2.get_time(t2_iterator)) { // Both TimeSeries have the current time point
                    if (add_values) {
                        merged_ts.add_time_point(ts1.get_time(t1_iterator),
                                                 ts1.get_value(t1_iterator) + ts2.get_value(t2_iterator));
                    }
                    else {
                        merged_ts.add_time_point(ts1.get_time(t1_iterator), ts1.get_value(t1_iterator));
                        log_warning("Both TimeSeries have values for t={}. The value of the first TimeSeries is used",
                                    ts1.get_time(t1_iterator));
                    }
                    t1_iterator += 1;
                    t2_iterator += 1;
                    if (t2_iterator >=
                        ts2.get_num_time_points()) { // Check if all values of second TimeSeries have been copied
                        t2_finished = true;
                        t2_iterator = ts2.get_num_time_points() - 1;
                    }
                }
                if (t1_iterator >=
                    ts1.get_num_time_points()) { // Check if all values of first TimeSeries have been copied
                    t1_finished = true;
                    t1_iterator = ts1.get_num_time_points() - 1;
                }
            }
            if (!t2_finished) {
                if (ts2.get_time(t2_iterator) < ts1.get_time(t1_iterator) ||
                    t1_finished) { // Current time point of second TimeSeries is smaller than current time point of first TimeSeries or first TimeSeries has already been copied entirely
                    merged_ts.add_time_point(ts2.get_time(t2_iterator), ts2.get_value(t2_iterator));
                    t2_iterator += 1;
                }
                else if (!t1_finished && ts2.get_time(t2_iterator) ==
                                             ts1.get_time(t1_iterator)) { // Both TimeSeries have the current time point
                    if (add_values) {
                        merged_ts.add_time_point(ts1.get_time(t1_iterator),
                                                 ts1.get_value(t1_iterator) + ts2.get_value(t2_iterator));
                    }
                    else {
                        merged_ts.add_time_point(ts1.get_time(t1_iterator), ts1.get_value(t1_iterator));
                        log_warning("Both TimeSeries have values for t={}. The value of the first TimeSeries is used",
                                    ts1.get_time(t1_iterator));
                    }
                    t1_iterator += 1;
                    t2_iterator += 1;
                    if (t1_iterator >=
                        ts1.get_num_time_points()) { // Check if all values of first TimeSeries have been copied
                        t1_finished = true;
                        t1_iterator = ts1.get_num_time_points() - 1;
                    }
                }
                if (t2_iterator >=
                    ts2.get_num_time_points()) { // Check if all values of second TimeSeries have been copied
                    t2_finished = true;
                    t2_iterator = ts2.get_num_time_points() - 1;
                }
            }
        }
    }
    return merged_ts;
}
} // namespace hybrid

} // namespace mio

#endif //MIO_HYBRID_ANALYZE_RESULT_H
