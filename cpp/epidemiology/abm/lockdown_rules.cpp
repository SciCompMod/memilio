/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#include "epidemiology/abm/lockdown_rules.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/time.h"


namespace epi
{

void set_home_office(TimePoint t_begin, double p, AbmMigrationParameters& params){
    auto damping1 = Eigen::VectorXd::Constant(1, p);
    params.get<WorkRatio>().add_damping(damping1, epi::SimulationTime(t_begin.days()));
}

void set_school_closure(TimePoint t_begin, double p, AbmMigrationParameters& params){
    auto damping1 = Eigen::VectorXd::Constant(1, p);
    params.get<SchoolRatio>().add_damping(damping1, epi::SimulationTime(t_begin.days()));
}

void close_social_events(TimePoint t_begin, double p, AbmMigrationParameters& params){
    auto damping1 = Eigen::VectorXd::Constant((size_t)AbmAgeGroup::Count, p);
    params.get<SocialEventRate>().add_damping(damping1, epi::SimulationTime(t_begin.days()));
}


} //namespace epi