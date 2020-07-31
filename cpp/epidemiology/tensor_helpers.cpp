#include "tensor_helpers.h"
#include <assert.h>
#include <algorithm>
#include <numeric>
#include <functional>

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
    //TODO: Maybe we need a version to process several flat indices in a batch?
    return unravel_index_given_prods(index, tensor_dimension_prods(dimensions));
}

std::vector<size_t> tensor_dimension_prods(std::vector<size_t> const& dimensions)
{
    // calculate [d_N, d_N*d_{N-1}, d_N*d_{N-1}*d_{N-2}, ...]
    std::vector<size_t> prods(dimensions.size());
    auto d   = dimensions.rbegin();
    prods[0] = *(d++);
    for (auto i = 1; i < prods.size() && d != dimensions.rend(); ++i) {
        prods[i] = prods[i - 1] * (*(d++));
    }
    return prods;
}

std::vector<size_t> unravel_index_given_prods(size_t const index, std::vector<size_t> const& dim_prods)
{
    assert(index >= 0 && index < dim_prods.back());

    std::vector<size_t> out(dim_prods.size(), 0);

    // recursively write out the remainder of the division by the dim_prods entries (in reverse order)
    int rem = index;
    int i   = 0;
    for (auto p = dim_prods.rbegin() + 1; p != dim_prods.rend(); ++p) {
        out[i++] = (int)rem / *p;
        rem      = rem % *p;
        if (rem == 0) {
            return out;
        }
    }
    out[dim_prods.size() - 1] = rem;

    return out;
}

namespace
{
    /**
     * @brief get_start_indices is an internal function used by `::get_slice_indices`
     *
     * It recursively gets all the start indices of all index ranges corresponding to one
     * slice of a hypothetical tensor of dimensions `dimensions`.
     *
     * @param current The current dyimension
     * @param target The target dimension of the desired slice
     * @param index The index into the target dimension of the desired slice
     * @param dimensions a vector of hypothetical dimension sizes of a the hypothetical tensor
     * @param starts A reference to the returned start indices
     * @param prod The product of all dimensions starting with dimension `current`.
     */
    void get_start_indices(size_t current, size_t target, size_t index, std::vector<size_t> const& dimensions,
                           std::vector<size_t>& starts, size_t& prod)
    {
        assert(dimensions[current] > 0);

        prod                    = prod / dimensions[current];
        size_t starts_init_size = starts.size();
        if (current < target) {
            if (current == 0) {
                for (size_t i = 0; i < dimensions[current]; ++i) {
                    starts.emplace_back(i * prod);
                }
                get_start_indices(current + 1, target, index, dimensions, starts, prod);
            }
            else {
                for (size_t j = 0; j < starts_init_size; ++j) {
                    for (size_t i = 1; i < dimensions[current]; ++i) {
                        starts.emplace_back(starts[j] + i * prod);
                    }
                }
                get_start_indices(current + 1, target, index, dimensions, starts, prod);
            }
        }
        else {
            if (target == 0) {
                starts.emplace_back(index * prod);
            }
            else {
                for (size_t j = 0; j < starts.size(); ++j) {
                    starts[j] += index * prod;
                }
            }
        }
    }

} // namespace

std::vector<size_t> get_slice_indices(size_t dimension, size_t index, std::vector<size_t> const& dimensions)
{
    // TODO: This should be reasonably fast if we are searching a long dimensions further to the left.
    // There might be more efficient ways to achieve this. Another (much simpler) option would be to
    // create a visitor that just iterates over all elements of a flat array and checks the unraveled index
    // against dimension and index using unravel_index.

    assert(dimension < dimensions.size());
    assert(dimensions[dimension] > 0);
    assert(index < dimensions[dimension]);

    size_t prod = std::accumulate(dimensions.begin(), dimensions.end(), 1, std::multiplies<size_t>());

    // recursively get all the start indices of the index ranges. prod is the size of the range
    std::vector<size_t> starts;
    starts.reserve(prod / dimensions[dimension]);
    get_start_indices(0, dimension, index, dimensions, starts, prod);

    // append the start indices by the ranges
    size_t starts_init_size = starts.size();
    for (size_t i = 0; i < starts_init_size; ++i) {
        std::vector<size_t> range(prod - 1);
        std::iota(range.begin(), range.end(), starts[i] + 1);
        starts.insert(starts.end(), range.begin(), range.end());
    }

    //sort and return
    std::sort(starts.begin(), starts.end());
    return starts;
}

} // namespace epi
