#include "common.h"
#include <assert.h>
namespace epi
{

size_t flatten_index(std::vector<size_t> const& indices, std::vector<size_t> const& dimensions)
{
    // calculate i_N + i_{N-1}*d_N + i_{N-2}*d_N*d_{N-1} + ... + i_1*d_N*d_{N-1}*...*d_2

    auto i     = indices.rbegin();
    size_t out = *(i++);

    auto d      = dimensions.rbegin();
    size_t prod = *(d++);

    while (i != indices.rend() && d != dimensions.rend()) {
        out += prod * (*(i++));
        prod *= *(d++);
    };

    return out;
}

std::vector<size_t> unravel_index(size_t const index, std::vector<size_t> const& dimensions)
{
    std::vector<size_t> out(dimensions.size());

    // pre-calculate [d_N, d_N*d_{N-1}, d_N*d_{N-1}*d_{N-2}, ...]
    //TODO: This is expensive. Maybe we need a version to process several flat indices in a batch?
    std::vector<size_t> prods(dimensions.size() - 1, 0);
    auto d   = dimensions.rbegin();
    prods[0] = *(d++);
    for (auto i = 1; i < prods.size() && d != dimensions.rend(); ++i) {
        prods[i] = prods[i - 1] * (*(d++));
    }

    assert(index >= 0 && index < prods.back()*dimensions[0]);

    // now write out the remainder of the division by the prod entries (in reverse order)
    int rem = index;
    int i   = 0;
    for (auto p = prods.rbegin(); p != prods.rend(); ++p) {
        out[i++] = (int)rem / *p;
        rem    = rem % *p;
        if (rem==0) {
            return out;
        }
    }
    out[dimensions.size() - 1] = rem;

    return out;
}

} // namespace epi
