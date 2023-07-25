/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Daniel Abele, Majid Abedi, Elisabeth Kluth
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
/** single include header for ABM */

#ifndef EPI_ABM_H
#define EPI_ABM_H

#include "abm/parameters.h"
#include "abm/simulation.h"
#include "abm/world.h"
#include "abm/person.h"
#include "abm/location.h"
#include "abm/location_type.h"
#include "memilio/utils/random_number_generator.h"
#include "abm/migration_rules.h"
#include "abm/testing_strategy.h"
#include "abm/infection.h"
#include "abm/infection_state.h"
#include "abm/virus_variant.h"
#include "abm/vaccine.h"
#include "abm/age.h"
#include "abm/household.h"
#include "abm/lockdown_rules.h"

#endif
