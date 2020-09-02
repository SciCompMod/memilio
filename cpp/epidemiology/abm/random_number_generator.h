#ifndef EPI_ABM_RANDOM_NUMBER_GENERATOR_H
#define EPI_ABM_RANDOM_NUMBER_GENERATOR_H

#include <random>
#include <functional>

namespace epi
{

/**
 * get a random number generator that is static and local to this thread.
 * @return a random number generator that is static and local to this thread.
 */
inline std::mt19937_64& thread_local_rng()
{
    static thread_local auto rng = std::mt19937_64(std::random_device()());
    return rng;
}

/**
 * adapter for a random number distribution.
 * Provides a static thread local instance of the distribution
 * and a replacable core generator function (e.g. for mocks during testing).
 * The parameters of the distribution are passed when the random number is generated
 * instead of when the object is constructed.
 * @tparam DistT a type that models the standard RandomNumberDistribution concept
 */
template<class DistT>
class DistributionAdapter
{
public:
    /**
     * The type returned by the distribution.
     */
    using ResultType = typename DistT::result_type;

    /**
     * The type that contains the parameters of the distribution.
     * The template parameter must be constructible from this type.
     */
    using ParamType = typename DistT::param_type;

    /**
     * The function that generates a random value from the distribution with the specified parameters.
     */
    using GeneratorFunction = std::function<ResultType(const ParamType& p)>;

    /**
     * the default generator function invokes an instance of the template parameter
     * with a static thread local RNG engine.
     */
    DistributionAdapter()
    {
        m_generator = [](auto&& params)
        {
            return DistT(params)(thread_local_rng());
        };
    }

    /**
     * get a random sample from the distribution.
     * accepts the same arguments as the constructors of the template parameter type.
     * example: 
     * std::uniform_int_distribution is constructed from two integers, so 
     * DistributionAdapter<std::uniform_int_distribution>::operator() accepts two integers as well.
     */
    template<class... T>
    ResultType operator()(T&&... params)
    {
        return m_generator(ParamType{std::forward<T>(params)...});
    }

    /**
     * get the generator function.
     */
    GeneratorFunction get_generator() const
    {
        return m_generator;
    }

    /**
     * set the generator function.
     */
    void set_generator(GeneratorFunction g)
    {
        m_generator = g;
    }

    /**
     * get a static thread local instance of this class.
     * Instance is default constructed on the first call.
     */
    static DistributionAdapter& get_instance()
    {
        static thread_local DistributionAdapter instance;
        return instance;
    }

private:
    GeneratorFunction m_generator;
};

/**
 * adapted std::exponential_distribution.
 * @see DistributionAdapter
 */
template<class Real>
using ExponentialDistribution = DistributionAdapter<std::exponential_distribution<Real>>;

/**
 * adapted std::uniform_int_distribution.
 * @see DistributionAdapter
 */
template<class Int>
using UniformIntDistribution = DistributionAdapter<std::uniform_int_distribution<Int>>;

/**
 * adapted uniform_real_distribution.
 * @see DistributionAdapter
 */
template<class Real>
using UniformDistribution = DistributionAdapter<std::uniform_real_distribution<Real>>;

} // namespace epi

#endif