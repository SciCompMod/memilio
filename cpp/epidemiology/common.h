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
 * @param indices a vector of indices of a tensor
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

} // namespace epi
