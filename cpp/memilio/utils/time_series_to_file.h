/*
* Copyright (C) 2023 Forschungszentrum Juelich - FZJ - IEK10, Germany
*
* Authors: Ralf Hannemann-Tamas
*
* Contact: Ralf Hannemann-tamas <r.hannemann-tamas@fz-juelich.de>
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
#ifndef EPI_TIME_SERIES_TO_FILE_H
#define EPI_TIME_SERIES_TO_FILE_H

#include "memilio/utils/time_series.h"
#include <fstream>


namespace mio {

/**
 * Print TimeSeries to file.
 * @tparam FP floating point time used by TimeSeries<FP>
 * @param time_series instance of time series to be printed
 * @param file_name name of the output file to be created
 */
template <class FP>
void time_series_to_file(const TimeSeries<FP>& time_series, const std::string& file_name) {
  std::ofstream output_file(file_name);
  output_file.precision(16);

  for(int i=0; i < time_series.get_num_time_points();++i) {
    for(int j=0; j < time_series.matrix().cols()-1; ++j) {
      output_file << time_series.matrix()(i,j) << " ";
      output_file << time_series.matrix()(i,time_series.matrix().cols()-1) << std::endl;
    }
  }

  output_file.close();
}

} // end namespace mio

#endif //EPI_TIME_SERIES_TO_FILE_H
