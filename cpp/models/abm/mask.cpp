#include "abm/mask_type.h"
#include "abm/parameters.h"
#include "abm/mask.h"
#include "abm/time.h"

namespace mio
{
namespace abm
{
Mask::Mask(MaskType type)
    : m_type(type)
    , m_time_used(TimeSpan(0))
{
}

void Mask::change_mask(MaskType new_mask_type)
{
    m_type      = new_mask_type;
    m_time_used = TimeSpan(0);
}

} // namespace abm
} // namespace mio