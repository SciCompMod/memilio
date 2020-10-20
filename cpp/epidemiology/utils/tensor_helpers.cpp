#include "epidemiology/utils/tensor_helpers.h"

#include <assert.h>
#include <algorithm>
#include <numeric>
#include <functional>

namespace epi
{

std::vector<size_t> unravel_index(size_t const index, std::vector<size_t> const& dimensions)
{
    //TODO: Maybe we need a version to process several flat indices in a batch?
    return unravel_index_given_prods(index, tensor_dimension_prods(dimensions));
}

std::vector<size_t> tensor_dimension_prods(std::vector<size_t> const& dimensions)
{
    // calculate [d_N, d_N*d_{N-1}, d_N*d_{N-1}*d_{N-2}, ...]
    std::vector<size_t> prods(dimensions.size());
    auto d   = dimensions.rbegin();
    prods[0] = *(d++);
    for (size_t i = 1; i < prods.size() && d != dimensions.rend(); ++i) {
        prods[i] = prods[i - 1] * (*(d++));
    }
    return prods;
}

std::vector<size_t> unravel_index_given_prods(size_t const index, std::vector<size_t> const& dim_prods)
{
    assert(index < dim_prods.back());

    std::vector<size_t> out(dim_prods.size(), 0);

    // recursively write out the remainder of the division by the dim_prods entries (in reverse order)
    size_t rem = index;
    size_t i   = 0;
    for (auto p = dim_prods.rbegin() + 1; p != dim_prods.rend(); ++p) {
        out[i++] = (size_t)rem / *p;
        rem      = rem % *p;
        if (rem == 0) {
            return out;
        }
    }
    out[dim_prods.size() - 1] = rem;

    return out;
}

} // namespace epi
