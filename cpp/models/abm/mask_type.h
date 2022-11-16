#ifndef EPI_ABM_MASK_TYPE_H
#define EPI_ABM_MASK_TYPE_H

#include <cstdint>

namespace mio
{
namespace abm
{

/**
 * type of a mask
 */
enum class MaskType : std::uint32_t
{
    Community = 0,
    Surgical,
    FFP2,

    Count
};
} // namespace abm
} // namespace mio

#endif