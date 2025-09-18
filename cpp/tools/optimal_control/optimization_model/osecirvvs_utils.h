#pragma once

#include <string>
#include <vector>

#include "models/ode_secirvvs/model.h"

using InfectionState = mio::osecirvvs::InfectionState;

static const std::vector<std::pair<std::string, InfectionState>> states_strings = {
    {"SusceptibleNaive", InfectionState::SusceptibleNaive},
    {"SusceptiblePartialImmunity", InfectionState::SusceptiblePartialImmunity},
    {"SusceptibleImprovedImmunity", InfectionState::SusceptibleImprovedImmunity},
    {"ExposedNaive", InfectionState::ExposedNaive},
    {"ExposedPartialImmunity", InfectionState::ExposedPartialImmunity},
    {"ExposedImprovedImmunity", InfectionState::ExposedImprovedImmunity},
    {"InfectedNoSymptomsNaive", InfectionState::InfectedNoSymptomsNaive},
    {"InfectedNoSymptomsPartialImmunity", InfectionState::InfectedNoSymptomsPartialImmunity},
    {"InfectedNoSymptomsImprovedImmunity", InfectionState::InfectedNoSymptomsImprovedImmunity},
    {"InfectedNoSymptomsNaiveConfirmed", InfectionState::InfectedNoSymptomsNaiveConfirmed},
    {"InfectedNoSymptomsPartialImmunityConfirmed", InfectionState::InfectedNoSymptomsPartialImmunityConfirmed},
    {"InfectedNoSymptomsImprovedImmunityConfirmed", InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed},
    {"InfectedSymptomsNaive", InfectionState::InfectedSymptomsNaive},
    {"InfectedSymptomsPartialImmunity", InfectionState::InfectedSymptomsPartialImmunity},
    {"InfectedSymptomsImprovedImmunity", InfectionState::InfectedSymptomsImprovedImmunity},
    {"InfectedSymptomsNaiveConfirmed", InfectionState::InfectedSymptomsNaiveConfirmed},
    {"InfectedSymptomsPartialImmunityConfirmed", InfectionState::InfectedSymptomsPartialImmunityConfirmed},
    {"InfectedSymptomsImprovedImmunityConfirmed", InfectionState::InfectedSymptomsImprovedImmunityConfirmed},
    {"InfectedSevereNaive", InfectionState::InfectedSevereNaive},
    {"InfectedSeverePartialImmunity", InfectionState::InfectedSeverePartialImmunity},
    {"InfectedSevereImprovedImmunity", InfectionState::InfectedSevereImprovedImmunity},
    {"InfectedCriticalNaive", InfectionState::InfectedCriticalNaive},
    {"InfectedCriticalPartialImmunity", InfectionState::InfectedCriticalPartialImmunity},
    {"InfectedCriticalImprovedImmunity", InfectionState::InfectedCriticalImprovedImmunity},
    {"DeadNaive", InfectionState::DeadNaive},
    {"DeadPartialImmunity", InfectionState::DeadPartialImmunity},
    {"DeadImprovedImmunity", InfectionState::DeadImprovedImmunity}};