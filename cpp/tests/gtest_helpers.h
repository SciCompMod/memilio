#ifndef GTEST_HELPERS_H
#define GTEST_HELPERS_H

#include <gmock/gmock.h>

//Action equivalent to assignment operator
//assignee and RHS value are free arguments
ACTION(Assign) { arg0 = arg1; }

//Action equivalent to assignment operator
//assignee is free argument, RHS value is fixed
ACTION_P(Assign, x) { arg0 = x; }

//Action equivalent to addition assignment operator
//assignee and RHS free arguments
ACTION(AddAssign) { arg0 += arg1; }

//Action equivalent to addition assignment operator
//assignee is free arguemtn, RHS value is fixed
ACTION_P(AddAssign, x) { arg0 += x; }

//versions that casts away const, necessary for fancy references like Eigen::Ref<Vector> to be assignable
//make sure arg0 is not originally defined const, otherwise undefined behaviour
ACTION_P(AddAssignUnsafe, x) { const_cast<arg0_type&>(arg0) += x; }
ACTION(AssignUnsafe) { const_cast<arg0_type&>(arg0) = arg1; }

#endif //GTEST_HELPERS_H