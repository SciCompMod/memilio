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
#ifndef EPI_DEREF_ITERATOR_H
#define EPI_DEREF_ITERATOR_H

#include "memilio/utils/transform_iterator.h"

namespace mio
{

namespace details
{
//helper functor that dereferences a pointer
struct DerefFunctor {
    template <class Ptr>
    auto& operator()(Ptr&& p) const
    {
        return *p;
    }
};

//base class for dereferencing iterator
template <class PtrIter>
using PointerDereferencingIteratorBase =
    decltype(make_transform_iterator(std::declval<PtrIter>(), std::declval<DerefFunctor>()));
} // namespace details

/**
 * iterator adaptor that makes a range of T* look like a range of T
 * e.g. std::vector<int*> v;
 * then *PointerDereferencingIterator<...>(v.begin()) == **v.begin()
 */
template <class PtrIter>
class PointerDereferencingIterator : public details::PointerDereferencingIteratorBase<PtrIter>
{
public:
    PointerDereferencingIterator(PtrIter ptr_iter)
        : details::PointerDereferencingIteratorBase<PtrIter>(ptr_iter)
    {
    }
};
} // namespace mio

#endif
