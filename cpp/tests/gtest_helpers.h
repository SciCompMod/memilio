#ifndef GTEST_HELPERS_H
#define GTEST_HELPERS_H

#include <gmock/gmock.h>

//Action equivalent to assignment operator
//first arg must be reference type
ACTION(Assign) { arg0 = arg1; }

//Action equivalent to addition assignment operator
//first arg must be reference type
ACTION(AddAssign) { arg0 += arg1; }

#endif //GTEST_HELPERS_H