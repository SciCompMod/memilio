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
#ifndef EPI_ABM_INFECTION_H
#define EPI_ABM_INFECTION_H

namespace mio
{
namespace abm
{

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
     * @return virus type of the infection
     */
    VirusType get_virus_type() const;
    
    
    InfectionState get_infection_state(const TimePoint t) const;
     
private:
    /**
     * determine viral load course and infection course
     */
    void draw_infection_course();
        //peak = gauss::get_instance();
    void determine_end_date();
    
    VirusType virus_type;
    TimePoint start_date;
    TimePoint end_date;
    double peak;
    double increase_slope;
    double decline_slope;
    double infectivity_parameter; // have to ask for distribution/parametrization of the infectivity
    //vector<pairs> Infection_course
    //TimePoint   state
    //start_date+2d         exposed
    //start_date+4d          carrier
    //...
    //start_date+xd           s/dead
    bool detected = false;
};

} // namespace abm
} // namespace mio

#endif
