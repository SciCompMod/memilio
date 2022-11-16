#include "abm/mask_type.h"
#include "abm/parameters.h"
#include "mask.h"

namespace mio
{
namespace abm
{
Mask::Mask(MaskType type)
    : m_type(type)
    , m_time_used(0)
{
}

double Mask::get_protection()
{
    if (m_type == MaskType::Community) {
        return 1.;
    }
    else if (m_type == MaskType::Surgical) {
        return 1.;
    }
    else if (m_type == MaskType::FFP2) {
        return 1.;
    }
    return 0;
}

void Mask::change_mask(MaskType new_mask_type)
{
    m_type = new_mask_type;
    m_time_used *= 0; // there's probably a better way
};

} // namespace abm
} // namespace mio