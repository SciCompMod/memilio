/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Kilian Volmer
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
#ifndef MIO_UTILS_BACK_INSERTER_SECOND_ELEMENT_H
#define MIO_UTILS_BACK_INSERTER_SECOND_ELEMENT_H

namespace mio
{
/**
* @brief Back inserter that ignores the first element of pairs given to it.
*
* Has the same member functions and objects as a std::back_insert_iterator.
* Template requires specialization for pair as it is not generally the same datatype as stored in Container.
*/
template <class Container, class Pair>
struct back_inserter_second_element {
    Container& container;
    back_inserter_second_element& operator*()
    {
        return *this;
    }
    back_inserter_second_element& operator++()
    {
        return *this;
    }
    back_inserter_second_element operator++(int)
    {
        return *this;
    }
    back_inserter_second_element operator=(const Pair& pair)
    {
        container.push_back(pair.second);
        return *this;
    }
};

} // namespace mio

#endif // MIO_UTILS_BACK_INSERTER_SECOND_ELEMENT_H
