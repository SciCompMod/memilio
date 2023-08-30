/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "memilio/utils/miompi.h"
#include "memilio/utils/span.h"

#include <cassert>
#include <functional>
#include <numeric>
#include <random>
#include <sstream>

namespace mio
{

/**
 * Models a uniform_random_bit_generator.
 * Keeps track of its seeds so they can be logged or set.
 * The generated sequence can be segmented into blocks that can help to reproduce 
 * simulations involving random numbers by reliably generating a specific part 
 * of the sequence or to assign each thread/process a different sequence.
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
        ++m_num_generated;
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

    /**
    * Seed this random number generator.
    * Starts a new random sequence at block 0.
    * @param seeds at least one seed, e.g., generated using std::random_device. More seeds increase quality.
    */
    void seed(const std::vector<unsigned int>& seeds)
    {
        assert(seeds.size() > 0);
        m_seeds = seeds;
        std::seed_seq sseq(m_seeds.begin(), m_seeds.end());
        m_rng.seed(sseq);
        m_num_generated = 0;
    }

    /**
    * Set the seeds in all MPI processes the same as in the root.
    */
    void synchronize_seeds()
    {
#ifdef MEMILIO_ENABLE_MPI
        int rank;
        MPI_Comm_rank(mpi::get_world(), &rank);
        int num_seeds;
        if (rank == 0) {
            num_seeds = int(m_seeds.size());
        }
        MPI_Bcast(&num_seeds, 1, MPI_INT, 0, mpi::get_world());
        if (rank != 0) {
            m_seeds.assign(num_seeds, 0);
        }
        MPI_Bcast(m_seeds.data(), num_seeds, MPI_UNSIGNED, 0, mpi::get_world());
#endif
    }

    /**
    * Get/Set the the size of blocks of the generated sequence.
    * This affects the number of samples skipped by forward_to_block().
    * Skipping samples is an O(n) operation, where n is the number of skipped samples, so choose the block size with care.
    * @{
    */
    void set_block_size(size_t block_size)
    {
        m_block_size = block_size;
    }
    size_t get_block_size() const
    {
        return m_block_size;
    }
    /**@}*/

    /**
    * Forward to block index i.
    * Skips numbers in the generated sequence up to the beginning of the block.
    * The next number generated will be element block_size * i of the random sequence.
    * This operation is O(n), where n is the number of skipped samples, so choose the block size with care.
    * @param i block index. May not be a block that is already passed or started.
    */
    void forward_to_block(size_t block_idx)
    {
        assert(block_idx * m_block_size >= m_num_generated &&
               "Can't forward to a previous block or one that is started.");
        auto num_remaining = m_block_size * block_idx - m_num_generated;
        m_rng.discard(num_remaining);
        m_num_generated += num_remaining;
    }

private:
    std::vector<unsigned int> m_seeds;
    std::mt19937_64 m_rng;
    size_t m_block_size = 1'000'000; ///< number of samples in a block, will be skipped by forward_to_block
    size_t m_num_generated = 0; ///< number of samples generated (or skipped) since the last seed
};

/**
 * get a random number generator that is static and local to this thread.
 * @return a random number generator that is static and local to this thread.
 */
RandomNumberGenerator& thread_local_rng();

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
    struct ParamType {
        using DistType = DistributionAdapter<DistT>;

        template <typename... Ps,
                  typename std::enable_if_t<std::is_constructible<typename DistT::param_type, Ps...>::value>* = nullptr>
        ParamType(Ps&&... ps)
            : params(std::forward<Ps>(ps)...)
        {
        }

        /**
        * get a static thread local instance of the contained Distribution class.
        * Calls DistributionAdapter::get_instance().
        */
        static DistributionAdapter& get_distribution_instance()
        {
            return DistType::get_instance();
        }

        typename DistT::param_type params;
    };
    /**
     * The function that generates a random value from the distribution with the specified parameters.
     */
    using GeneratorFunction = std::function<ResultType(const typename DistT::param_type& p)>;

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
        return m_generator(typename DistT::param_type{std::forward<T>(params)...});
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
