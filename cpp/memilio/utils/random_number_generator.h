/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef EPI_ABM_RANDOM_NUMBER_GENERATOR_H
#define EPI_ABM_RANDOM_NUMBER_GENERATOR_H

#include "memilio/utils/logging.h"
#include "memilio/utils/span.h"

#ifdef MEMILIO_HAS_MPI
#include "mpi.h"
#endif

#include <cassert>
#include <functional>
#include <numeric>
#include <random>
#include <sstream>

namespace mio
{

/**
 * models a uniform_random_bit_generator.
 * keeps track of its seeds so they can be logged or set.
 * @see thread_local_rng for a static instance.
 */
class RandomNumberGenerator
{
public:
    using result_type = std::mt19937_64::result_type;

    static constexpr result_type min()
    {
        return std::mt19937_64::min();
    }
    static constexpr result_type max()
    {
        return std::mt19937_64::max();
    }
    result_type operator()()
    {
        return m_rng();
    }

    static std::vector<unsigned int> generate_seeds()
    {
        std::random_device rd;
        return {rd(), rd(), rd(), rd(), rd(), rd()};
    }

    RandomNumberGenerator()
        : m_seeds(generate_seeds())
    {
        std::seed_seq sseq(m_seeds.begin(), m_seeds.end());
        m_rng.seed(sseq);
    }
    std::vector<unsigned int> get_seeds() const
    {
        return m_seeds;
    }
    void seed(const std::vector<unsigned int>& seeds)
    {
        m_seeds = seeds;
        std::seed_seq sseq(m_seeds.begin(), m_seeds.end());
        m_rng.seed(sseq);
    }

private:
    std::vector<unsigned int> m_seeds;
    std::mt19937_64 m_rng;
};

/**
 * models a uniform_random_bit_generator, using cached numbers instead of generating on demand.
 * @see thread_local_rng for a static instance when MPI is enabled.
 */
class CachedRNG : public mio::RandomNumberGenerator
{
public:
    using result_type = mio::RandomNumberGenerator::result_type;

    CachedRNG(const size_t cache_size, const std::vector<result_type>& predrawn_samples)
        : RandomNumberGenerator()
        , m_counter(0)
        , m_index(0)
        , m_cache(predrawn_samples)
    {
        assert(cache_size > 0);
        m_cache.resize(cache_size);
        // fill remaining cache using RandomNumberGenerator::operator()
        for (size_t i = predrawn_samples.size(); i < m_cache.size(); i++) {
            m_cache[i] = (*this)();
        }
    };

    ~CachedRNG()
    {
        mio::log_info("CachedRNG used {} out of {} cached numbers, looped cache {} time(s).",
                      m_index + m_counter * (long long int)m_cache.size(), m_cache.size(), m_counter);
    }

    result_type operator()()
    {
        if (m_index >= m_cache.size()) {
            m_counter++;
            m_index = 0;
            mio::log_warning("CachedRNG has too few cached numbers! The cache will be looped on subsequent calls. Increase the cache size to avoid degraded RNG.");
        }
        return m_cache[m_index++];
    }
    size_t cache_size() const
    {
        return m_cache.size();
    }

protected:
    size_t m_counter, m_index;
    std::vector<result_type> m_cache;
};

#ifndef MIO_GLOBAL_RANDOM_NUMBER_GENERATOR
    #define MIO_GLOBAL_RANDOM_NUMBER_GENERATOR
    #ifdef MEMILIO_HAS_MPI
        /**
         * @brief get a random number generator using cached numbers that is static and local to this thread.
         * 
         * The first call to this function can be used to set the cache size, and optionally predrawn samples.
         * All subsequent calls ignore these parameters, and only return the thread_local CachedRNG instance.
         * The cache is reused after running out of numbers, logging a warning. To avoid degraded random numbers,
         * increase the cache size if such a warning occurs.
         * 
         * @param cache_size number of random numbers to be drawn and stored.  
         * @param samples 
         * @return a random number generator that is static and local to this thread.
         */
        inline CachedRNG& thread_local_rng(const size_t cache_size=100000, const std::vector<mio::CachedRNG::result_type>& samples={})
        {
            static thread_local auto rng = CachedRNG(cache_size, samples);
            return rng;
        }
    #else
        /**
         * get a random number generator that is static and local to this thread.
         * @return a random number generator that is static and local to this thread.
         */
        inline RandomNumberGenerator& thread_local_rng()
        {
            static thread_local auto rng = RandomNumberGenerator();
            return rng;
        }
    #endif
#endif

#ifdef MEMILIO_HAS_MPI
/**
 * @brief Draws local_cache_size random numbers for each rank on comm from the root rank, and distributes them.
 *
 * This function must be called after MPI_Init and before any call to mio::thread_local_rng.
 * Draws numbers using a RandomNumberGenerator instance on root, then scatters them across all ranks on comm.
 * Expects the default definition of thread_local_rng with MEMILIO_HAS_MPI enabled.
 *
 * @param local_cache_size The size of the cache on each rank. Must be the same across all ranks.
 * @param seeds Vector of seeds for the RandomNumberGenerator.
 * @param root The rank that draws and distibutes all random numbers.
 * @param comm MPI Communicator where the numbers should be distributed on.
 */
void initialize_rng_cache(const size_t local_cache_size, const std::vector<unsigned int>& seeds={}, const int root=0,
                          const MPI_Comm comm=MPI_COMM_WORLD)
{
    int my_rank, my_size;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &my_size);

    std::vector<CachedRNG::result_type> buffer, rng_cache(local_cache_size);
    RandomNumberGenerator rng;
    if (seeds.size() != 0) {
        rng.seed(seeds);
    }
    if (my_rank == root) {
        buffer.resize(local_cache_size * my_size);
        for (auto& rn : buffer)
            rn = rng();
    }
    MPI_Scatter(buffer.data(), local_cache_size, MPI_UINT64_T, rng_cache.data(), local_cache_size, MPI_UINT64_T, 0,
                comm);
    mio::thread_local_rng(local_cache_size, rng_cache);
}
#endif

inline void log_rng_seeds(const RandomNumberGenerator& rng, LogLevel level)
{
    const auto& seeds = rng.get_seeds();
    std::stringstream ss;
    bool first = true;
    for (auto& s : seeds) {
        if (!first) {
            ss << ", ";
        }
        first = false;
        ss << s;
    }
    log(level, "Using RNG with seeds: {0}.", ss.str());
}

inline void log_thread_local_rng_seeds(LogLevel level)
{
    log_rng_seeds(thread_local_rng(), level);
}

/**
 * adapter for a random number distribution.
 * Provides a static thread local instance of the distribution
 * and a replacable core generator function (e.g. for mocks during testing).
 * The parameters of the distribution are passed when the random number is generated
 * instead of when the object is constructed.
 * @tparam DistT a type that models the standard RandomNumberDistribution concept
 */
template <class DistT>
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
        m_generator = [](auto&& params) {
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
    template <class... T>
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
 * select a random integer in [0, n) with weights [w_0, ..., w_(n-1)]
 * the probability to pick i is w_i/S where S is the sum of all weights.
 * similar to std::discrete_distribution but does not allocate.
 * Models the standard RandomNumberDistribution concept.
 */
template <class Int>
class DiscreteDistributionInPlace
{
public:
    /**
     * the type returned by the distribution.
     */
    using result_type = Int;

    /**
     * stores the parameters of the distribution (i.e. the weights).
     */
    class param_type
    {
    public:
        using distribution = DiscreteDistributionInPlace;

        param_type() = default;

        param_type(Span<double> weights)
            : m_weights(weights)
        {
        }

        Span<double> weights() const
        {
            return m_weights;
        }

    private:
        Span<double> m_weights;
    };

    /**
     * default distribution has no weights.
     * always returns 0.
     */
    DiscreteDistributionInPlace() = default;

    /**
     * distribution with specified weights.
     */
    DiscreteDistributionInPlace(Span<double> weights)
        : m_params(weights)
    {
    }

    /**
     * distribution with specified params.
     */
    DiscreteDistributionInPlace(param_type params)
        : m_params(params)
    {
    }

    /**
     * reset internal state.
     * does nothing, but required by the concept.
     */
    void reset()
    {
    }

    /**
     * get the parameters of the distribution.
     */
    param_type param()
    {
        return m_params;
    }

    /**
     * set the parameters of the distribution.
     */
    void param(param_type p)
    {
        m_params = p;
    }

    /**
     * get the weights.
     */
    Span<double> weights()
    {
        return m_params.weights();
    }

    /**
     * draw a random number from the distribution.
     * @param rng object of a type that that models UniformRandomBitGenerator concept.
     */
    template <class RNG>
    result_type operator()(RNG& rng)
    {
        return (*this)(rng, m_params);
    }

    /**
     * draw a random number from the distribution with the specified parameters.
     * @param rng object of a type that that models UniformRandomBitGenerator concept.
     * @param p parameters of the dstribution.
     */
    template <class RNG>
    result_type operator()(RNG& rng, param_type p)
    {
        auto weights = p.weights();
        if (weights.size() <= 1) {
            return 0;
        }
        auto sum = std::accumulate(weights.begin(), weights.end(), 0.0);
        auto u =
            std::uniform_real_distribution<double>()(rng, std::uniform_real_distribution<double>::param_type{0.0, sum});
        auto intermediate_sum = 0.0;
        for (size_t i = 0; i < weights.size(); ++i) {
            intermediate_sum += weights.get_ptr()[i];
            if (u < intermediate_sum) {
                return i;
            }
        }
        assert(false && "this should never happen.");
        return result_type(-1);
    }

private:
    param_type m_params;
};

/**
 * adapted discrete distribution
 * @see DistributionAdapter
 */
template <class Int>
using DiscreteDistribution = DistributionAdapter<DiscreteDistributionInPlace<Int>>;

/**
 * adapted std::exponential_distribution.
 * @see DistributionAdapter
 */
template <class Real>
using ExponentialDistribution = DistributionAdapter<std::exponential_distribution<Real>>;

/**
 * adapted std::uniform_int_distribution.
 * @see DistributionAdapter
 */
template <class Int>
using UniformIntDistribution = DistributionAdapter<std::uniform_int_distribution<Int>>;

/**
 * adapted uniform_real_distribution.
 * @see DistributionAdapter
 */
template <class Real>
using UniformDistribution = DistributionAdapter<std::uniform_real_distribution<Real>>;

} // namespace mio

#endif
