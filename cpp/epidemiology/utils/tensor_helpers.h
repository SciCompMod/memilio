#ifndef TENSOR_HELPERS_H
#define TENSOR_HELPERS_H

#include <assert.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <functional>

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
template <class Container>
void get_start_indices(size_t current, size_t target, size_t index, Container const& dimensions,
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

namespace epi
{

/**
 * @brief flatten_index takes a set of indices into a mutlidemsional array and calculates the flat index
 *
 * Given indices (i,j,k,...) of a tensor with dimensions (n,m,l,...), flatten_index calculates
 * the index of the corresponding element if the elements are sorted sequentially in a
 * row major fashion (that is right indices are incremented before left indices)
 *
 * @param indices a vector of indices of a hypothetical tensor
 * @param dimensions a vector of the dimension sizes of each dimension
 * @return the corresponding flat index
 */
template <class IndexContainer = std::initializer_list<size_t>, class DimContainer>
size_t flatten_index(IndexContainer const& indices, DimContainer const& dimensions)
{

#ifndef NDEBUG
    assert(indices.size() == dimensions.size());
    size_t j = 0;
    for (auto idx : indices) {
        assert(idx < dimensions[j++]);
    }
#endif

    // calculate i_N + i_{N-1}*d_N + i_{N-2}*d_N*d_{N-1} + ... + i_1*d_N*d_{N-1}*...*d_2

    auto i     = std::rbegin(indices);
    size_t out = *(i++);

    auto d      = dimensions.rbegin();
    size_t prod = *(d++);

    while (i != std::rend(indices) && d != dimensions.rend()) {
        out += prod * (*(i++));
        prod *= *(d++);
    };

    return out;
}

/**
 * @brief unravel_index takes a flat index of a tensor and calculates the index of each dimension
 *
 * If the input index is an index into a sequentially stored multidimensional array of dimensions
 * (n,m,l,...), unravel_index calculcates the corresponding index for each dimension. It is assumed
 * that the array is stored in row major fashion (that is right indices are incremented before
 * left indices)
 *
 * @param index the flat index to be converted
 * @param dimensions a vector of the dimension sizes of each dimension
 * @return a vector of indices
 */
std::vector<size_t> unravel_index(size_t const index, std::vector<size_t> const& dimensions);

/**
 * @brief tensor_dimension_prods given the dimensions (d_0,...,d_N) of a hypothetical tensor,
 * `::tensor_dimension_prods` returns the vector of dimension products
 *
 *    [d_N, d_N*d_{N-1}, d_N*d_{N-1}*d_{N-2}, ..., d_N*...*d_0].
 *
 * This vector can be used e.g. in `::unravel_index_given_prods`.
 *
 * @param dimensions a vector of the dimension sizes of each dimension
 * @return a vector of dimension size products
 */
std::vector<size_t> tensor_dimension_prods(std::vector<size_t> const& dimensions);

/**
 * @brief unravel_index_given_prods is the same as `::unravel_index` with a different input for better performance
 *
 * `::unravel_index` calculates the index of each dimension given a flat index into a hypothetical tensor
 * of given dimensions (d_0, ..., d_N). `::unravel_index` computes the vector of dimension products
 *
 *   dim_prods = [d_N, d_N*d_{N-1}, d_N*d_{N-1}*d_{N-2}, ..., d_N*...*d_0]
 *
 * for the dimensions (d_0, ...., d_N). Instead of the vector of dimensions, `::unravel_index_given_prods`
 * takes the precalculated vector dim_prods as an input, so that it does not have to be recalculated with each
 * call, but can be cached somewhere else. The vector can be calculated using `::tensor_dimension_prods`.
 *
 * @param indices  a vector of indices of a hypothetical tensor
 * @param dim_prods a vector containing the dimension products
 * @return the corresponding flat index
 */
std::vector<size_t> unravel_index_given_prods(size_t const index, std::vector<size_t> const& dim_prods);

/**
 * @brief get_slice_indices returns all the flat indices corresponding to an index into one dimension of a
 * hypothetical tensor. It is assumed that the tensor is stored in row-major fashion.
 *
 * Example: If dimensions=(3, 3), we have a hypothetical 2rd-order tesor of size 3x3. Stored in row-major
 * fashion the flat indices are
 *
 *    0 : (0, 0)
 *    1 : (0, 1)
 *    2 : (0, 2)
 *    3 : (1, 0)
 *    4 : (1, 1)
 *    5 : (1, 2)
 *    6 : (2, 0)
 *    7 : (2, 1)
 *    8 : (2, 2)
 *
 * `get_slice_indices(0, 1, {3,3})` will return all indices corresponding to the first dimension with index==1,
 * that is the indices of the second row (3,4,5)
 *
 * `get_slice_indices(1, 0, {3,3})` will return the indices (0,3,6) corresponding to the second dimension
 * with index==0, i.e. the first column.
 *
 * @param dimension the dimension of the slice
 * @param index the index of the slice
 * @param dimensions a vector of dimension sizes of a hypothetical tensor
 * @return the flat indices of all elements in the desired slice
 */
template <typename Container>
std::vector<size_t> get_slice_indices(size_t dimension, size_t index, Container const& dimensions)
{
    // TODO: This should be reasonably fast if we are searching a long dimensions further to the left.
    // There might be more efficient ways to achieve this. Another (much simpler) option would be to
    // create a visitor that just iterates over all elements of a flat array and checks the unraveled index
    // against dimension and index using unravel_index.
    // an even better option would be to use boost::multi_array and the implemented views

    assert(dimension < dimensions.size());
    assert(dimensions[dimension] > 0);
    assert(index < dimensions[dimension]);

    size_t prod = std::accumulate(dimensions.begin(), dimensions.end(), size_t(1), std::multiplies<size_t>());

    // recursively get all the start indices of the index ranges. prod is the size of the range
    std::vector<size_t> starts;
    starts.reserve(prod / dimensions[dimension]);
    get_start_indices(0, dimension, index, dimensions, starts, prod);
    std::sort(starts.begin(), starts.end());

    // append the start indices by the ranges
    for (auto it = starts.begin(); it != starts.end();) {
        starts.insert(it + 1, prod - 1, 0);
        std::iota(it + 1, it + prod, *it + 1);
        it += prod;
    }

    return starts;
}

} // namespace epi

#endif // TENSOR_HELPERS_H
