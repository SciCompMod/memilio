/* 
* Copyright (C) 2020-2024 MEmilio
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

// Core ABM components
#include "abm/parameters.h"
#include "abm/simulation.h"
#include "abm/world.h"
#include "abm/person.h"
#include "abm/location.h"
#include "abm/location_type.h"
#include "abm/infection.h"
#include "abm/infection_state.h"
#include "abm/virus_variant.h"
#include "abm/time.h"

// Movement and migration
#include "abm/migration_rules.h"
#include "abm/trip_list.h"
#include "abm/movement_data.h"

// Testing and interventions
#include "abm/testing_strategy.h"
#include "abm/test_type.h"
#include "abm/vaccine.h"
#include "abm/mask.h"
#include "abm/mask_type.h"
#include "abm/lockdown_rules.h"

// Social structures
#include "abm/household.h"

// Analysis and utilities
#include "abm/analyze_result.h"
#include "abm/random_events.h"
#include "abm/config.h"

// Mathematical utilities
#include "memilio/math/interpolation.h"
#include "memilio/utils/random_number_generator.h"
#endif
