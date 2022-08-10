/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: David Kerkmann
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
#ifndef EPI_ABM_LOCATION_H
#define EPI_ABM_LOCATION_H

namespace mio
{
namespace abm
{

class Infectivity{

public:
    /**
     * create an infectivity curve
     * @param virus_type specified virus type
     */
    Infectivity(VirusType virus_type);

    /**
     * get infectivity at a given time
     * @param t time point of the querry
     * @return infectivity at given time point
     */
    double get_infectivity(const TimePoint t) const;
    
private:
    void draw_infection_course();
    
    double slope; // have to ask for distribution/parametrization of the infectivity
}

class Infection{
    
public:
    /**
     * create an infection for a single person.
     * @param start_date starting date of the infection
     * @param virus_type virus type of the infection
     */
    Infection(TimePoint start_date, VirusType virus_type);
    
    /**
     * get viral load at a given time
     * @param t time point of the querry
     * @return viral load at given time point
     */
    double get_viral_load(const TimePoint t) const;
    
    /**
     * get infectivity at a given time
     * @param t time point of the querry
     * @return infectivity at given time point
     */
    double get_infectivity(const TimePoint t) const;
    
    /**
     * get virus type
     */
    const VirusType get_virus_type() const;
    
    
private:
    void draw_infection_course();
    void determine_end_date();
    
    VirusType virus_type;
    TimePoint start_date;
    TimePoint end_date;
    double peak;
    double increase_slope;
    double decline_slope;
}

} // namespace abm
} // namespace mio

#endif
