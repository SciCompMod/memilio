/* 
* Copyright (C) 2020-2025 MEmilio
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

#ifndef MIO_RANDOM_NUMBER_GENERATOR_H
#define MIO_RANDOM_NUMBER_GENERATOR_H

#include "memilio/io/default_serialize.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/span.h"
#include "memilio/utils/type_safe.h"

MSVC_WARNING_DISABLE_PUSH(4127) //conditional expression is constant
GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wexpansion-to-defined") //Random123 handles the portability of this warning internally
#include "Random123/array.h"
#include "Random123/threefry.h"
GCC_CLANG_DIAGNOSTIC(pop)
MSVC_WARNING_POP()

#include <cassert>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <type_traits>

namespace mio
{

/**
* Base class for counter based random number generator.
*
* All (pseudo) random number generators (RNG) consist of some state, 
* a function `state = advance(state)`, and a function `sample = generate(state)`.
* They produce a sequence of random samples by advancing the state and from the new state 
* generate the actual sample. 
*
* example in pseudo code:
* state = initial_state(seed)
* sample1 = generate(state)
* state = advance(state)
* sample2 = generate(state)
* ...
*
* In most normal RNGs like mersenne twister the advance function needs to
* be complicated in order to perform the necessary mixing of bits. The state
* needs to be relatively large to contain sufficient numbers of random bits. The generate function
* on the other hand is very simple, often just returning the state (whole or in parts) without further
* computation. 
*
* In counter based generators (cRNG), the state and advance function are very simple and 
* all the mixing of bits is in the generate function. The state is split into a key and a counter.
* The generate function is an encryption or hash function. Like a hash function, it produces a
* pseudo-random value from the input. The key is used by the generate function the same way as an encryption key.
* The key is randomly seeded in the beginning so that a different sequence of samples is generated on each run.
* Then the key doesn't change anymore. The advance function only increments the counter. The generate function
* produces completely different output even for sequential numbers.
*
* Because their state is simple, cRNG are well suited for parallel applications. To create n independent subsequences from
* a total sequence of N samples you only need to create n counters, where counter i starts at i * (N / n).
* Normal RNGs need special algorithms to efficiently generate subsequences, if it is possible at all.
* The counter is of minimal size, it only needs to be big enough to fit the number of samples generated.
* Often the subsequence index i is already available in some other form, e.g., the thread index or the agent index in an 
* agent based model, so only a small amount of extra storage is needed for the subsequence counter. 
* The key is shared between all subsequences. Modern CPU architectures also are very efficient at executing the hash and
* encryption functions that are used as the generate function, increasing performance of the generator. 
* 
* The length of the total sequence and the subsequences and the number of subsequences can be adjusted as needed
* by assigning different numbers of bits to the key, the subsequence index and the subsequence counter.
* The counter only needs to store the number of samples generated. A counter of c bits supports a 
* sequence of 2^c samples or 2^n subsequences of 2^(c - n) samples, in which case you can split the 
* counter into a subsequence index of n bits and a subsequence counter of (c - n) bits where each
* of the subsequence counters starts at 0.
* Generating a samples of k bits requires a key of at least k bits for sufficient randomness.
* Example:
* * A 64 bit counter (uint64_t) and a 64 bit key produce 2^64 samples each with 64 bits. You need to store
* one counter and one key.
* * A subsequence index of 32 bits (uint32_t), subsequence counter with 32 bits and a 64 bit key
* produce 2^32 subsequences of 2^32 samples each with 64 bits per sample. You need to store 2^32 subsequence indices, 
* 2^32 subsequence counters and one key, but the counters are completely independent and thread safe. The subsequence index
* and corresponding subsequence counter can also be stored together in one 64bit counter, e.g., the subsequence index is the 
* high bits and subsequence counter is the low bits, see `rng_totalsequence_index()`.
*
* Also see https://github.com/DEShawResearch/random123 for more information on cRNGs and the specific cRNG
* we use.
* 
* Classes deriving from this base class need to supply the key and counter by implementing
* the functions
* * result_type get_key() const
* * result_type get_counter() const
* * void increment_counter()
*
* This class satisfies the standard UniformRandomBitGenerator concept.
*/
template <class Derived>
class RandomNumberGeneratorBase
{
public:
    using result_type = uint64_t;

    /**
    * Minimum value generated by this generator.
    * Counterbased generators allow the whole range supported by the result_type.
    */
    static constexpr result_type min()
    {
        return std::numeric_limits<result_type>::min();
    }

    /**
    * Maximum value generated by this generator.
    * Counterbased generators allow the whole range supported by the result_type.
    */
    static constexpr result_type max()
    {
        return std::numeric_limits<result_type>::max();
    }

    /**
    * Generate the next value in the random sequence.
    * Key and counter are supplied by the Derived class by implementing
    * the functions get_key(), get_counter(), and increment_counter().
    */
    result_type operator()();
};

namespace details
{
/**
* Convert a Random123 array type (rng counters and keys) to uint64_t.
*/
inline uint64_t to_uint64(r123array2x32 tf_array)
{
    uint64_t i;
    std::memcpy(&i, tf_array.data(), sizeof(uint64_t));
    return i;
}

/**
* Convert a uint64_t to a Random123 array type (rng counters and keys).
*/
inline r123array2x32 to_r123_array(uint64_t i)
{
    threefry2x32_ctr_t c;
    std::memcpy(c.data(), &i, sizeof(uint64_t));
    return c;
}
} // namespace details

/**
* A key type for counter based random number generators.
* @tparam an unsigned integer type that determines the size of the key, i.e., the number of different sequences.
*/
template <class T>
struct MEMILIO_ENABLE_EBO Key : TypeSafe<T, Key<T>>, OperatorComparison<Key<T>> {
    static_assert(std::is_unsigned<T>::value, "Underlying Integer type must be unsigned.");
    using TypeSafe<T, Key<T>>::TypeSafe;
};
static_assert(sizeof(Key<uint32_t>) == sizeof(uint32_t), "Empty Base Optimization isn't working.");

/**
* A counter type for counter based random number generators.
* @tparam an unsigned integer type that determines the size of the counter, i.e., the length of the random sequence.
*/
template <class T>
struct MEMILIO_ENABLE_EBO Counter : TypeSafe<T, Counter<T>>,
                                    OperatorComparison<Counter<T>>,
                                    OperatorAdditionSubtraction<Counter<T>> {
    static_assert(std::is_unsigned<T>::value, "Underlying Integer type must be unsigned.");
    using TypeSafe<T, Counter<T>>::TypeSafe;
};
static_assert(sizeof(Counter<uint32_t>) == sizeof(uint32_t), "Empty Base Optimization isn't working.");

template <class Derived>
auto RandomNumberGeneratorBase<Derived>::operator()() -> result_type
{
    //generate a random sample using the Random123 library.
    //Use another threefryNxR algorithm if larger or more random samples are needed than 64 bit.
    auto self = static_cast<Derived*>(this);
    auto c    = static_cast<uint64_t>(self->get_counter().get());
    auto k    = static_cast<uint64_t>(self->get_key().get());
    auto r    = details::to_uint64(threefry2x32(details::to_r123_array(k), details::to_r123_array(c)));
    self->increment_counter();
    return r;
}

/**
* Seed a counter based random number generator key.
* @tparam A type that satisfies standard SeedSeq.
* @param seed_seq A seed sequence, e.g. initialized using std::random_device.
* @return A seeded key.
*/
template <class SeedSeq>
Key<uint64_t> seed_rng_key(SeedSeq& seed_seq)
{
    auto tf_key = threefry2x32_key_t::seed(seed_seq);
    return Key<uint64_t>(details::to_uint64(tf_key));
}

/**
* Get the counter in the total sequence for a counter in a given subsequence.
* A total sequence counter of C bits supports 2^N subsequences of 2^(C - N) samples. Then the counter 
* can be split into a subsequence index of N bits and a subsequence counter of S = (C - N) bits.
* The length of the subsequences is determined by the type of the subsequence counter and the 
* requested total sequence counter.
* The subsequence index does not need to be exactly N bits, it can be larger. E.g., a 
* subsequence counter of 16 bits and a subsequence index of 64 bits can be used for
* a combined counter of 64 bits since there is not 48 bit type.
* @tparam UIntC The counter type of the total sequence with C bits.
* @tparam UIntN An unsigned integer type with at least N bits for the subsequence index.
* @tparam CounterS A counter type with S = (C - N) bits for the subsequence counter.
* @param subsequence_idx The index of the subsequence. Must be less than 2^N.
* @param counter The counter in the subsequence.
* @param return The counter in the total sequence.
*/
template <class UIntC, class UIntN, class CounterS>
Counter<UIntC> rng_totalsequence_counter(UIntN subsequence_idx, CounterS counter)
{
    //use UIntC for variables because it's the biggest integer type in this function
    static const UIntC BITS_PER_BYTE = 8;
    static const UIntC C_BITS        = sizeof(UIntC) * BITS_PER_BYTE;
    static const UIntC S_BITS        = sizeof(CounterS) * BITS_PER_BYTE;
    static const UIntC N_BITS        = C_BITS - S_BITS;

    static_assert(S_BITS < C_BITS, "Subsequence counter must be smaller than total sequence counter.");
    static_assert(N_BITS <= C_BITS, "Subsequence index must not be bigger than total sequence counter.");
    static_assert(N_BITS <= sizeof(UIntN) * BITS_PER_BYTE, "Subsequence index must be at least N bits");

    assert(UIntC(subsequence_idx) <= (UIntC(1) << N_BITS) &&
           "Subsequence index is too large."); //(1 << N) is the same as (2^N)

    //N high bits: subsequence idx
    //S low bits: subsequence counter
    //=> C = N + S bits: total sequence counter
    //example:
    //subsequence index uint32_t(181) = 0x000000B5
    //subsequence counter uint32_t(41309) = 0x0000A15D
    //total sequence counter = 0x000000B50000A15D
    const auto i = static_cast<UIntC>(subsequence_idx);
    const auto s = static_cast<UIntC>(counter.get());
    const auto c = (i << S_BITS) + s; //shift subsequence index to the high bits, add subsequence counter into low bits
    return Counter<UIntC>{c};
}

/**
* Get the subsequence counter from the total sequence counter.
* A total sequence counter of C bits supports 2^N subsequences of 2^(C - N) samples. Then the counter 
* can be split into a subsequence index of N bits and a subsequence counter of S = (C - N) bits.
* The length of the subsequences is determined by the type of the subsequence counter.
* @tparam UIntS An unsigned integer type of S bits for the subsequence counter.
* @tparam CounterC A counter type of C bits, where C > S.
* @param counter The total sequence counter.
* @return The counter in the subsequence.
*/
template <class UIntS, class CounterC>
Counter<UIntS> rng_subsequence_counter(CounterC counter)
{
    using UIntC                = typename CounterC::ValueType;
    static const UIntC C_BYTES = sizeof(UIntC);
    static const UIntC S_BYTES = sizeof(UIntS);

    static_assert(S_BYTES < C_BYTES, "Subsequence counter must be smaller than total sequence counter.");

    //the subsequence counter is in the lower bits of total sequence counter
    //see rng_totalsequence_counter above
    //so we just need to truncate the counter to get the subsequence counter
    return Counter<UIntS>(static_cast<UIntS>(counter.get()));
}

/**
* General purpose counter based random number generator.
* Stores its own key and counter.
* @see RandomNumberGeneratorBase.
*/
class RandomNumberGenerator : public RandomNumberGeneratorBase<RandomNumberGenerator>
{
public:
    RandomNumberGenerator()
        : m_counter(0)
    {
        seed(generate_seeds());
    }

    Key<uint64_t> get_key() const
    {
        return m_key;
    }
    Counter<uint64_t> get_counter() const
    {
        return m_counter;
    }
    void set_counter(Counter<uint64_t> counter)
    {
        m_counter = counter;
    }
    void increment_counter()
    {
        ++m_counter;
    }
    static std::vector<uint32_t> generate_seeds()
    {
        std::random_device rd;
        return {rd(), rd(), rd(), rd(), rd(), rd()};
    }

    void seed(const std::vector<uint32_t>& seeds)
    {
        m_seeds = seeds;
        std::seed_seq seed_seq(m_seeds.begin(), m_seeds.end());
        m_key = seed_rng_key(seed_seq);
    }

    const std::vector<uint32_t> get_seeds() const
    {
        return m_seeds;
    }

    /**
    * Set the seeds in all MPI processes the same as in the root.
    */
    void synchronize()
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
        if (rank != 0) {
            seed(m_seeds);
        }
#endif
    }

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("RandomNumberGenerator").add("key", m_key).add("counter", m_counter).add("seeds", m_seeds);
    }

private:
    Key<uint64_t> m_key;
    Counter<uint64_t> m_counter;
    std::vector<uint32_t> m_seeds;
};

/**
 * get a random number generator that is static and local to this thread.
 * @return a RandomNumberGenerator that is static and local to this thread.
 * @note Not to be used anymore, only used by ParameterDistribution.
 */
RandomNumberGenerator& thread_local_rng();

/**
* Log the seeds used by the RandomNumberGenerator at the specified LogLevel.
*/
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

/**
* Log the seeds used by the RandomNumberGenerator from thread_local_rng() at the specified LogLevel.
*/
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

private:
    /**
     * the default generator function invokes an instance of the template parameter
     * with a static thread local RNG engine.
     * Constructors are private, use get_instance to get the current version.
     */
    DistributionAdapter()                                      = default;
    DistributionAdapter(const DistributionAdapter&)            = default;
    DistributionAdapter& operator=(const DistributionAdapter&) = default;
    DistributionAdapter(DistributionAdapter&&)                 = default;
    DistributionAdapter& operator=(DistributionAdapter&&)      = default;

public:
    /**
     * get a random sample from the distribution.
     * accepts the same arguments as the constructors of the template parameter type.
     * example: 
     * std::uniform_int_distribution is constructed from two integers, so 
     * DistributionAdapter<std::uniform_int_distribution>::operator() accepts two integers as well.
     */
    template <class RNG, class... T>
    ResultType operator()(RNG& rng, T&&... params)
    {
        if (m_generator) {
            //unlikely outside of tests
            return m_generator(typename DistT::param_type{std::forward<T>(params)...});
        }
        else {
            return DistT(std::forward<T>(params)...)(rng);
        }
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
     * get a static instance of this class.
     * Instance is default constructed on the first call.
     * The generator function of this instance can be replaced
     * for mocking during tests.
     */
    static DistributionAdapter& get_instance()
    {
        static DistributionAdapter instance;
        return instance;
    }

private:
    GeneratorFunction m_generator;
};

/**
 * select a random integer in [0, n) with weights [w_0, ..., w_(n-1)]
 * the probability to pick i is w_i/S where S is the sum of all weights.
 * Similar to std::discrete_distribution but does not allocate, instead
 * expects a Span of weights.
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

template <class IOContext, class UniformDistributionParams,
          class Real              = typename UniformDistributionParams::DistType::ResultType,
          std::enable_if_t<std::is_same_v<UniformDistributionParams, typename UniformDistribution<Real>::ParamType>,
                           void*> = nullptr>
void serialize_internal(IOContext& io, const UniformDistributionParams& p)
{
    auto obj = io.create_object("UniformDistributionParams");
    obj.add_element("a", p.params.a());
    obj.add_element("b", p.params.b());
}

template <class IOContext, class UniformDistributionParams,
          class Real              = typename UniformDistributionParams::DistType::ResultType,
          std::enable_if_t<std::is_same_v<UniformDistributionParams, typename UniformDistribution<Real>::ParamType>,
                           void*> = nullptr>
IOResult<UniformDistributionParams> deserialize_internal(IOContext& io, Tag<UniformDistributionParams>)
{
    auto obj = io.expect_object("UniformDistributionParams");
    auto a   = obj.expect_element("a", Tag<Real>{});
    auto b   = obj.expect_element("b", Tag<Real>{});
    return apply(
        io,
        [](auto&& a_, auto&& b_) {
            return UniformDistributionParams{a_, b_};
        },
        a, b);
}

/**
 * adapted poisson_distribution.
 * @see DistributionAdapter
 */
template <class Int>
using PoissonDistribution = DistributionAdapter<std::poisson_distribution<Int>>;

/**
 * adapted lognormal_distribution.
 * @see DistributionAdapter
 */
template <class Real>
using LogNormalDistribution = DistributionAdapter<std::lognormal_distribution<Real>>;

} // namespace mio

#endif
