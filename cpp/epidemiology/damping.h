#ifndef DAMPING_H
#define DAMPING_H

#include <vector>

namespace epi
{

/**
 * This defined a damping factor for a
 * mitigation strategy for one point in time.
 */
class Damping
{
public:
    double day;
    double factor;

    Damping(double day_in, double factor_in);
};

class Dampings
{
public:
    /**
     * @brief default constructor, contains constant damping of 1 everywhere.
     */
    Dampings();

    /**
     * @brief activates or deactivates smoothing between two distinct discrete dampings
     */
    void set_smoothing(bool smoothing);


    /**
     * @brief Adds a damping to the current model
     * @param d The damping, which is a factor and day from which the mitigation acts
     */
    void add(const Damping& d);

    /**
     * @brief Returns the damping factor
     *
     * @param[in] day Current day
     * @param[in] smooth smoothing dampings by cosine curve to avoid stepping function
     */
    double get_factor(double day) const;

private:
    std::vector<Damping> m_dampings;
    bool m_smoothing;
};

} // namespace epi

#endif // DAMPING_H
