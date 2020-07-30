#pragma once

#include <vector>

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
size_t flatten_index(std::vector<size_t> const& indices, std::vector<size_t> const& dimensions);

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

} // namespace epi
