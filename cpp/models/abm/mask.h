#include "abm/mask_type.h"
#include "abm/time.h"

namespace mio
{
namespace abm
{
class Mask
{
public:
    Mask(MaskType type);

    MaskType get_type()
    {
        return m_type;
    }

    TimeSpan get_time_used()
    {
        return m_time_used;
    }

    void increase_time_used(TimeSpan dt)
    {
        m_time_used += dt;
    }

    double get_protection(); // dependend on the type

    void change_mask(MaskType new_mask_type); // changes the type of the mask and sets time_used to 0

private:
    MaskType m_type;
    TimeSpan m_time_used;
};
} // namespace abm
} // namespace mio