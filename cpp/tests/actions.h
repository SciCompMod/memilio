/*
* Copyright (C) 2020-2026 MEmilio
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
#ifndef EPI_TESTS_ACTIONS_H
#define EPI_TESTS_ACTIONS_H

#include <gmock/gmock.h>

//Action equivalent to assignment operator
//assignee and RHS value are free arguments
ACTION(Assign)
{
    arg0 = arg1;
}

//Action equivalent to assignment operator
//assignee is free argument, RHS value is fixed
ACTION_P(Assign, x)
{
    arg0 = x;
}

//Action equivalent to addition assignment operator
//assignee and RHS free arguments
ACTION(AddAssign)
{
    arg0 += arg1;
}

//Action equivalent to addition assignment operator
//assignee is free arguemtn, RHS value is fixed
ACTION_P(AddAssign, x)
{
    arg0 += x;
}

//versions that casts away const, necessary for fancy references like Eigen::Ref<Vector> to be assignable
//make sure arg0 is not originally defined const, otherwise undefined behaviour
ACTION_P(AddAssignUnsafe, x)
{
    const_cast<std::decay_t<arg0_type>&>(arg0) += x;
}
ACTION(AssignUnsafe)
{
    const_cast<std::decay_t<arg0_type>&>(arg0) = arg1;
}

#endif //EPI_TESTS_ACTIONS_H
