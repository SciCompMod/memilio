#include "models/abm/activity_type.h"
#include "models/abm/infection_state.h"

#include <vector>

#ifndef JULIA_TEST_PKG_H
#define JULIA_TEST_PKG_H

#ifdef __cplusplus
extern "C" {
#endif

struct MyStruct {
    int a;
    double b;
};

struct InfectionTuple {
    int time_point; // time point of aquiring the infection state in seconds since the start of the simulation
    mio::abm::InfectionState infection_state;
};

struct PersonData {
    int id, age_group;
    uint32_t activity_type;
};

void printint(int n);
float mean_of_two(int x, int y);

int struct_test(MyStruct s);

bool hospital_model(PersonData p, const InfectionTuple* infection_history, size_t history_len, int time_of_day);

#ifdef __cplusplus
}
#endif

#endif
