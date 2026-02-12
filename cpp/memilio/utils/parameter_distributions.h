/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Martin J. Kuehn, Daniel Abele, Julia Bicker
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
#ifndef PARAMETER_DISTRIBUTIONS_H
#define PARAMETER_DISTRIBUTIONS_H

#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/visitor.h"
#include "memilio/utils/random_number_generator.h"
#include "models/abm/personal_rng.h"
#include "memilio/io/io.h"

#include <limits>
#include <memory>
#include <vector>
#include <random>

namespace mio
{

/**
 * @brief This is a visitor class to visit all Parameter Distribution objects
 *
 * More information to the visitor pattern is here: https://en.wikipedia.org/wiki/Visitor_pattern
 */
using ParameterDistributionVisitor =
    Visitor<class ParameterDistributionNormal, class ParameterDistributionUniform, class ParameterDistributionLogNormal,
            class ParameterDistributionExponential, class ParameterDistributionConstant>;
using ConstParameterDistributionVisitor =
    ConstVisitor<class ParameterDistributionNormal, class ParameterDistributionUniform,
                 class ParameterDistributionLogNormal, class ParameterDistributionExponential,
                 class ParameterDistributionConstant>;

template <class Derived>
using VisitableParameterDistribution =
    Visitable<Derived, class ParameterDistribution, ParameterDistributionVisitor, ConstParameterDistributionVisitor>;

namespace details
{
template <class IOObj>
struct SerializationVisitor : ConstParameterDistributionVisitor {
    SerializationVisitor(IOObj& o)
        : obj(o)
    {
    }
    virtual void visit(const ParameterDistributionNormal& normal_dist) final;
    virtual void visit(const ParameterDistributionUniform& uniform_dist) final;
    virtual void visit(const ParameterDistributionLogNormal& lognormal_dist) final;
    virtual void visit(const ParameterDistributionExponential& lognormal_dist) final;
    virtual void visit(const ParameterDistributionConstant& lognormal_dist) final;
    IOObj& obj;
};
} // namespace details

/*
 * Parameter Distribution class representing a generic distribution and contains predefined samples
 */
class ParameterDistribution
{
public:
    ParameterDistribution()
    {
    }

    virtual ~ParameterDistribution() = default;

    void add_predefined_sample(ScalarType sample)
    {
        m_predefined_samples.push_back(sample);
    }

    void remove_predefined_samples()
    {
        m_predefined_samples.resize(0);
    }

    const std::vector<ScalarType>& get_predefined_samples() const
    {
        return m_predefined_samples;
    }

    /*
     * @brief returns a value for the given parameter distribution
     * in case some predefined samples are set, these values are taken
     * first, in case the vector of predefined values is empty, a 'real'
     * random sample is taken
     */
    template <class RNG>
    ScalarType get_sample(RNG& rng)
    {
        if (m_predefined_samples.size() > 0) {
            ScalarType rnumb = m_predefined_samples[0];
            m_predefined_samples.erase(m_predefined_samples.begin());
            return rnumb;
        }
        else {
            return get_rand_sample(rng);
        }
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj     = io.create_object("ParameterDistribution");
        auto visitor = details::SerializationVisitor<decltype(obj)>{obj};
        this->accept(visitor);
    }

    virtual bool smaller_impl(const ParameterDistribution& other) const = 0;

    /**
     * @brief Returns the distribution parameters as vector.
     */
    virtual std::vector<ScalarType> params() const = 0;

    virtual ScalarType get_rand_sample(RandomNumberGenerator& rng)          = 0;
    virtual ScalarType get_rand_sample(abm::PersonalRandomNumberGenerator&) = 0;

    virtual ParameterDistribution* clone() const = 0;

    /**
     * @brief This function implements the visit interface of the visitor pattern
     *
     * It can be used for any ways of working with the class to dispatch
     * the class type. More information here: https://en.wikipedia.org/wiki/Visitor_pattern
     */
    virtual void accept(ParameterDistributionVisitor& visitor)            = 0;
    virtual void accept(ConstParameterDistributionVisitor& visitor) const = 0;

protected:
    std::vector<ScalarType>
        m_predefined_samples; // if these values are set; no real sample will occur but these values will be taken
};

/*
 * Child class of Parameter Distribution class which additionally contains
 * the mean value and the standard deviation for a normal distribution
 */
class ParameterDistributionNormal : public VisitableParameterDistribution<ParameterDistributionNormal>
{
public:
    ParameterDistributionNormal()
        : VisitableParameterDistribution<ParameterDistributionNormal>()
        , m_mean(0)
        , m_standard_dev(1)
        , m_distribution(0, 1)
    {
    }

    ParameterDistributionNormal(ScalarType mean, ScalarType standard_dev)
        : VisitableParameterDistribution<ParameterDistributionNormal>()
        , m_mean(mean)
        , m_standard_dev(standard_dev)
        , m_distribution(mean, standard_dev)
    {
        m_mean         = mean;
        m_standard_dev = standard_dev;
    }

    ParameterDistributionNormal(ScalarType lower_bound, ScalarType upper_bound, ScalarType mean, ScalarType quantile)
        : VisitableParameterDistribution<ParameterDistributionNormal>()
        , m_mean(mean)
        , m_upper_bound(upper_bound)
        , m_lower_bound(lower_bound)
        , m_quantile(quantile)
    {
        // if upper and lower bound are given, the standard deviation is calculated such that [lower_bound, upper_bound] represent the 0.995 quartile]
        m_standard_dev = upper_bound; // set as to high and adapt then
        adapt_standard_dev(m_standard_dev, upper_bound, lower_bound);
        m_distribution = mio::NormalDistribution<ScalarType>::ParamType(m_mean, m_standard_dev);
    }

    ParameterDistributionNormal(ScalarType lower_bound, ScalarType upper_bound, ScalarType mean,
                                ScalarType standard_dev, ScalarType quantile)
        : VisitableParameterDistribution<ParameterDistributionNormal>()
        , m_mean(mean)
        , m_standard_dev(standard_dev)
        , m_upper_bound(upper_bound)
        , m_lower_bound(lower_bound)
        , m_quantile(quantile)
    {
        check_quantiles(m_mean, m_standard_dev);
        m_distribution = mio::NormalDistribution<ScalarType>::ParamType(m_mean, m_standard_dev);
    }

    void set_mean(ScalarType mean)
    {
        m_mean = mean;
    }

    bool check_quantiles()
    {
        return check_quantiles(m_mean, m_standard_dev);
    }

    /*
     * @brief verification that at least 99% of the density
     * function lie in the interval defined by the boundaries
     */
    bool check_quantiles(ScalarType& mean, ScalarType& standard_dev)
    {
        bool changed = false;
        if (adapt_mean(mean)) {
            changed = true;
            log_warning("Mean adapted to lie within [lowerbound,upperbound].");
        }

        if (adapt_standard_dev(standard_dev, m_upper_bound, m_lower_bound)) {
            changed = true;
            log_warning("Standard deviation reduced to fit 99% of the distribution within [lowerbound,upperbound].");
        }

        return changed;
    }

    bool adapt_mean(ScalarType& mean)
    {
        bool changed = false;
        if (mean < m_lower_bound || mean > m_upper_bound) {
            mean    = 0.5 * (m_upper_bound - m_lower_bound);
            changed = true;
        }
        return changed;
    }

    // ensure that 0.99 % of the distribution are within lower bound and upper bound
    bool adapt_standard_dev(ScalarType& standard_dev, ScalarType upper_bound, ScalarType lower_bound)
    {
        bool changed = false;
        if (m_mean + standard_dev * m_quantile > upper_bound) {
            standard_dev = (upper_bound - m_mean) / m_quantile;
            changed      = true;
        }
        if (m_mean - standard_dev * m_quantile < lower_bound) {
            standard_dev = (m_mean - lower_bound) / m_quantile;
            changed      = true;
        }

        return changed;
    }

    void set_standard_dev(ScalarType standard_dev)
    {
        m_standard_dev = standard_dev;
    }

    ScalarType get_mean() const
    {
        return m_mean;
    }

    ScalarType get_standard_dev() const
    {
        return m_standard_dev;
    }

    void log_stddev_changes(bool log_stddev_change)
    {
        m_log_stddev_change = log_stddev_change;
    }

    void set_lower_bound(ScalarType lower_bound)
    {
        m_lower_bound = lower_bound;
    }

    void set_upper_bound(ScalarType upper_bound)
    {
        m_upper_bound = upper_bound;
    }

    ScalarType get_lower_bound() const
    {
        return m_lower_bound;
    }

    ScalarType get_upper_bound() const
    {
        return m_upper_bound;
    }

    bool smaller_impl(const ParameterDistribution& other) const override
    {
        const ParameterDistributionNormal* PDNormal = dynamic_cast<const ParameterDistributionNormal*>(&other);
        if (PDNormal) {
            return m_mean < PDNormal->m_mean;
        }
        else {
            return false;
        }
    }

    std::vector<ScalarType> params() const override
    {
        return {m_mean, m_standard_dev};
    }

    /*
     * @brief gets a sample of a normally distributed variable
     * before sampling, it is verified that at least 99% of the
     * density function lie in the interval defined by the boundaries
     * otherwise the normal distribution is adapted
     */
    template <class RNG>
    ScalarType sample(RNG& rng)
    {
        //If ub = lb, sampling can only be succesful if mean = lb and dev = 0.
        //But this degenerate normal distribution is not allowed by the c++ standard.
        if (m_upper_bound == m_lower_bound) {
            return m_lower_bound;
        }

        if (check_quantiles(m_mean, m_standard_dev) || m_distribution.params.mean() != m_mean ||
            m_distribution.params.stddev() != m_standard_dev) {
            m_distribution = NormalDistribution<ScalarType>::ParamType{m_mean, m_standard_dev};
        }

        int i            = 0;
        int retries      = 10;
        ScalarType rnumb = m_distribution.get_distribution_instance()(thread_local_rng(), m_distribution.params);
        while ((rnumb > m_upper_bound || rnumb < m_lower_bound) && i < retries) {
            rnumb = m_distribution.get_distribution_instance()(rng, m_distribution.params);
            i++;
            if (i == retries) {
                log_warning("Not successfully sampled within [min,max].");
                if (rnumb > m_upper_bound) {
                    rnumb = m_upper_bound;
                }
                else {
                    rnumb = m_lower_bound;
                }
            }
        }
        return rnumb;
    }

    ScalarType get_rand_sample(RandomNumberGenerator& rng) override
    {
        return sample(rng);
    }

    ScalarType get_rand_sample(abm::PersonalRandomNumberGenerator& rng) override
    {
        return sample(rng);
    }

    template <class IOObject>
    void serialize_elements(IOObject& obj) const
    {
        obj.add_element("Mean", m_mean);
        obj.add_element("StandardDev", m_standard_dev);
        obj.add_element("LowerBound", m_lower_bound);
        obj.add_element("UpperBound", m_upper_bound);
        obj.add_element("Quantile", m_quantile);
        obj.add_list("PredefinedSamples", m_predefined_samples.begin(), m_predefined_samples.end());
    }

    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("ParameterDistributionNormal");
        serialize_elements(obj);
    }

    template <class IOContext, class IOObject>
    static IOResult<ParameterDistributionNormal> deserialize_elements(IOContext& io, IOObject& obj)
    {
        auto m      = obj.expect_element("Mean", Tag<ScalarType>{});
        auto s      = obj.expect_element("StandardDev", Tag<ScalarType>{});
        auto lb     = obj.expect_element("LowerBound", Tag<ScalarType>{});
        auto ub     = obj.expect_element("UpperBound", Tag<ScalarType>{});
        auto qu     = obj.expect_element("Quantile", Tag<ScalarType>{});
        auto predef = obj.expect_list("PredefinedSamples", Tag<ScalarType>{});
        auto p      = apply(
            io,
            [](auto&& lb_, auto&& ub_, auto&& m_, auto&& s_, auto&& qu_, auto&& predef_) {
                auto distr = ParameterDistributionNormal(lb_, ub_, m_, s_, qu_);
                for (auto&& e : predef_) {
                    distr.add_predefined_sample(e);
                }
                return distr;
            },
            lb, ub, m, s, qu, predef);
        if (p) {
            return success(p.value());
        }
        else {
            return p.as_failure();
        }
    }

    template <class IOContext>
    static IOResult<ParameterDistributionNormal> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("ParameterDistributionNormal");
        return deserialize_elements(io, obj);
    }

    ParameterDistribution* clone() const override
    {
        return new ParameterDistributionNormal(*this);
    }

private:
    ScalarType m_mean; // the mean value of the normal distribution
    ScalarType m_standard_dev; // the standard deviation of the normal distribution
    ScalarType m_upper_bound = std::numeric_limits<
        ScalarType>::max(); // upper bound and lower bound can be given to the constructor instead of stddev
    ScalarType m_lower_bound = std::numeric_limits<ScalarType>::min();
    ScalarType m_quantile    = 2.5758; // default is 0.995 quartile
    NormalDistribution<ScalarType>::ParamType m_distribution;
    bool m_log_stddev_change = true;
};

template <class IOObj>
void details::SerializationVisitor<IOObj>::visit(const ParameterDistributionNormal& normal_dist)
{
    obj.add_element("Type", std::string("Normal"));
    normal_dist.serialize_elements(obj);
}

/*
 * Child class of Parameter Distribution class which represents an uniform distribution 
 */
class ParameterDistributionUniform : public VisitableParameterDistribution<ParameterDistributionUniform>
{
public:
    ParameterDistributionUniform(ScalarType lower_bound, ScalarType upper_bound)
        : VisitableParameterDistribution<ParameterDistributionUniform>()
        , m_upper_bound(upper_bound)
        , m_lower_bound(lower_bound)
        , m_distribution(lower_bound, upper_bound)
    {
    }

    bool smaller_impl(const ParameterDistribution& other) const override
    {
        const ParameterDistributionUniform* PDUniform = dynamic_cast<const ParameterDistributionUniform*>(&other);
        if (PDUniform) {
            return (m_lower_bound < PDUniform->m_lower_bound) && (m_upper_bound < PDUniform->m_upper_bound);
        }
        else {
            return false;
        }
    }

    std::vector<ScalarType> params() const override
    {
        return {m_lower_bound, m_upper_bound};
    }

    ScalarType get_lower_bound() const
    {
        return m_lower_bound;
    }

    ScalarType get_upper_bound() const
    {
        return m_upper_bound;
    }

    /*
     * @brief gets a sample of a uniformly distributed variable
     */
    template <class RNG>
    ScalarType sample(RNG& rng)
    {
        if (m_distribution.params.b() != m_upper_bound || m_distribution.params.a() != m_lower_bound) {
            m_distribution = UniformDistribution<ScalarType>::ParamType{m_lower_bound, m_upper_bound};
        }

        return m_distribution.get_distribution_instance()(rng, m_distribution.params);
    }

    ScalarType get_rand_sample(RandomNumberGenerator& rng) override
    {
        return sample(rng);
    }

    ScalarType get_rand_sample(abm::PersonalRandomNumberGenerator& rng) override
    {
        return sample(rng);
    }

    ParameterDistribution* clone() const override
    {
        return new ParameterDistributionUniform(*this);
    }

    template <class IOObject>
    void serialize_elements(IOObject& obj) const
    {
        obj.add_element("LowerBound", m_lower_bound);
        obj.add_element("UpperBound", m_upper_bound);
        obj.add_list("PredefinedSamples", m_predefined_samples.begin(), m_predefined_samples.end());
    }

    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("ParameterDistributionUniform");
        serialize_elements(obj);
    }

    template <class IOContext, class IOObject>
    static IOResult<ParameterDistributionUniform> deserialize_elements(IOContext& io, IOObject& obj)
    {
        auto lb     = obj.expect_element("LowerBound", Tag<ScalarType>{});
        auto ub     = obj.expect_element("UpperBound", Tag<ScalarType>{});
        auto predef = obj.expect_list("PredefinedSamples", Tag<ScalarType>{});
        auto p      = apply(
            io,
            [](auto&& lb_, auto&& ub_, auto&& predef_) {
                auto distr = ParameterDistributionUniform(lb_, ub_);
                for (auto&& e : predef_) {
                    distr.add_predefined_sample(e);
                }
                return distr;
            },
            lb, ub, predef);
        if (p) {
            return success(p.value());
        }
        else {
            return p.as_failure();
        }
    }

    template <class IOContext>
    static IOResult<ParameterDistributionUniform> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("ParameterDistributionUniform");
        return deserialize_elements(io, obj);
    }

private:
    ScalarType m_upper_bound;
    ScalarType m_lower_bound;
    UniformDistribution<ScalarType>::ParamType m_distribution;
};

template <class IOObj>
void details::SerializationVisitor<IOObj>::visit(const ParameterDistributionUniform& uniform_dist)
{
    obj.add_element("Type", std::string("Uniform"));
    uniform_dist.serialize_elements(obj);
}

/*
 * Child class of Parameter Distribution class which represents an lognormal distribution 
 */
class ParameterDistributionLogNormal : public VisitableParameterDistribution<ParameterDistributionLogNormal>
{
public:
    ParameterDistributionLogNormal(ScalarType log_mean, ScalarType log_stddev)
        : VisitableParameterDistribution<ParameterDistributionLogNormal>()
        , m_log_mean(log_mean)
        , m_log_stddev(log_stddev)
        , m_distribution(log_mean, log_stddev)
    {
    }

    bool smaller_impl(const ParameterDistribution& other) const override
    {
        const ParameterDistributionLogNormal* PDLogNorm = dynamic_cast<const ParameterDistributionLogNormal*>(&other);
        if (PDLogNorm) {
            return m_log_mean < PDLogNorm->m_log_mean;
        }
        else {
            return false;
        }
    }

    std::vector<ScalarType> params() const override
    {
        return {m_log_mean, m_log_stddev};
    }

    ScalarType get_log_mean() const
    {
        return m_log_mean;
    }

    ScalarType get_log_stddev() const
    {
        return m_log_stddev;
    }

    /*
     * @brief gets a sample of a lognormally distributed variable
     */
    template <class RNG>
    ScalarType sample(RNG& rng)
    {
        if (m_distribution.params.m() != m_log_mean || m_distribution.params.s() != m_log_stddev) {
            m_distribution = LogNormalDistribution<ScalarType>::ParamType{m_log_mean, m_log_stddev};
        }

        return m_distribution.get_distribution_instance()(rng, m_distribution.params);
    }

    ScalarType get_rand_sample(RandomNumberGenerator& rng) override
    {
        return sample(rng);
    }

    ScalarType get_rand_sample(abm::PersonalRandomNumberGenerator& rng) override
    {
        return sample(rng);
    }

    ParameterDistribution* clone() const override
    {
        return new ParameterDistributionLogNormal(*this);
    }

    template <class IOObject>
    void serialize_elements(IOObject& obj) const
    {
        obj.add_element("LogMean", m_log_mean);
        obj.add_element("LogStddev", m_log_stddev);
        obj.add_list("PredefinedSamples", m_predefined_samples.begin(), m_predefined_samples.end());
    }

    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("ParameterDistributionLogNormal");
        serialize_elements(obj);
    }

    template <class IOContext, class IOObject>
    static IOResult<ParameterDistributionLogNormal> deserialize_elements(IOContext& io, IOObject& obj)
    {
        auto lm     = obj.expect_element("LogMean", Tag<ScalarType>{});
        auto ls     = obj.expect_element("LogStddev", Tag<ScalarType>{});
        auto predef = obj.expect_list("PredefinedSamples", Tag<ScalarType>{});
        auto p      = apply(
            io,
            [](auto&& lm_, auto&& ls_, auto&& predef_) {
                auto distr = ParameterDistributionLogNormal(lm_, ls_);
                for (auto&& e : predef_) {
                    distr.add_predefined_sample(e);
                }
                return distr;
            },
            lm, ls, predef);
        if (p) {
            return success(p.value());
        }
        else {
            return p.as_failure();
        }
    }

    template <class IOContext>
    static IOResult<ParameterDistributionLogNormal> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("ParameterDistributionLogNormal");
        return deserialize_elements(io, obj);
    }

private:
    ScalarType m_log_mean;
    ScalarType m_log_stddev;
    LogNormalDistribution<ScalarType>::ParamType m_distribution;
};

template <class IOObj>
void details::SerializationVisitor<IOObj>::visit(const ParameterDistributionLogNormal& uniform_dist)
{
    obj.add_element("Type", std::string("LogNormal"));
    uniform_dist.serialize_elements(obj);
}

/*
 * Child class of Parameter Distribution class which represents an exponential distribution 
 */
class ParameterDistributionExponential : public VisitableParameterDistribution<ParameterDistributionExponential>
{
public:
    ParameterDistributionExponential(ScalarType rate)
        : VisitableParameterDistribution<ParameterDistributionExponential>()
        , m_rate(rate)
        , m_distribution(rate)
    {
    }

    bool smaller_impl(const ParameterDistribution& other) const override
    {
        const ParameterDistributionExponential* PDExponential =
            dynamic_cast<const ParameterDistributionExponential*>(&other);
        if (PDExponential) {
            return m_rate < PDExponential->m_rate;
        }
        else {
            return false;
        }
    }

    std::vector<ScalarType> params() const override
    {
        return {m_rate};
    }

    ScalarType get_rate() const
    {
        return m_rate;
    }

    /*
     * @brief gets a sample of a exponentially distributed variable
     */
    template <class RNG>
    ScalarType sample(RNG& rng)
    {
        if (m_distribution.params.lambda() != m_rate) {
            m_distribution = ExponentialDistribution<ScalarType>::ParamType{m_rate};
        }

        return m_distribution.get_distribution_instance()(rng, m_distribution.params);
    }

    ScalarType get_rand_sample(RandomNumberGenerator& rng) override
    {
        return sample(rng);
    }

    ScalarType get_rand_sample(abm::PersonalRandomNumberGenerator& rng) override
    {
        return sample(rng);
    }

    ParameterDistribution* clone() const override
    {
        return new ParameterDistributionExponential(*this);
    }

    template <class IOObject>
    void serialize_elements(IOObject& obj) const
    {
        obj.add_element("Rate", m_rate);
        obj.add_list("PredefinedSamples", m_predefined_samples.begin(), m_predefined_samples.end());
    }

    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("ParameterDistributionExponential");
        serialize_elements(obj);
    }

    template <class IOContext, class IOObject>
    static IOResult<ParameterDistributionExponential> deserialize_elements(IOContext& io, IOObject& obj)
    {
        auto r      = obj.expect_element("Rate", Tag<ScalarType>{});
        auto predef = obj.expect_list("PredefinedSamples", Tag<ScalarType>{});
        auto p      = apply(
            io,
            [](auto&& r_, auto&& predef_) {
                auto distr = ParameterDistributionExponential(r_);
                for (auto&& e : predef_) {
                    distr.add_predefined_sample(e);
                }
                return distr;
            },
            r, predef);
        if (p) {
            return success(p.value());
        }
        else {
            return p.as_failure();
        }
    }

    template <class IOContext>
    static IOResult<ParameterDistributionExponential> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("ParameterDistributionExponential");
        return deserialize_elements(io, obj);
    }

private:
    ScalarType m_rate;
    ExponentialDistribution<ScalarType>::ParamType m_distribution;
};

template <class IOObj>
void details::SerializationVisitor<IOObj>::visit(const ParameterDistributionExponential& uniform_dist)
{
    obj.add_element("Type", std::string("Exponential"));
    uniform_dist.serialize_elements(obj);
}

/*
 * Child class of Parameter Distribution class which represents a constant distribution/value 
 */
class ParameterDistributionConstant : public VisitableParameterDistribution<ParameterDistributionConstant>
{
public:
    ParameterDistributionConstant(ScalarType constant)
        : VisitableParameterDistribution<ParameterDistributionConstant>()
        , m_constant(constant)
    {
    }

    bool smaller_impl(const ParameterDistribution& other) const override
    {
        const ParameterDistributionConstant* PDConstantl = dynamic_cast<const ParameterDistributionConstant*>(&other);
        if (PDConstantl) {
            return m_constant < PDConstantl->m_constant;
        }
        else {
            return false;
        }
    }

    std::vector<ScalarType> params() const override
    {
        return {m_constant};
    }

    ScalarType get_constant() const
    {
        return m_constant;
    }

    /*
     * @brief gets a constant
     */
    template <class RNG>
    ScalarType sample(RNG& /*rng*/)
    {
        return m_constant;
    }

    ScalarType get_rand_sample(RandomNumberGenerator& rng) override
    {
        return sample(rng);
    }

    ScalarType get_rand_sample(abm::PersonalRandomNumberGenerator& rng) override
    {
        return sample(rng);
    }

    ParameterDistribution* clone() const override
    {
        return new ParameterDistributionConstant(*this);
    }

    template <class IOObject>
    void serialize_elements(IOObject& obj) const
    {
        obj.add_element("Constant", m_constant);
        obj.add_list("PredefinedSamples", m_predefined_samples.begin(), m_predefined_samples.end());
    }

    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("ParameterDistributionConstant");
        serialize_elements(obj);
    }

    template <class IOContext, class IOObject>
    static IOResult<ParameterDistributionConstant> deserialize_elements(IOContext& io, IOObject& obj)
    {
        auto c      = obj.expect_element("Constant", Tag<ScalarType>{});
        auto predef = obj.expect_list("PredefinedSamples", Tag<ScalarType>{});
        auto p      = apply(
            io,
            [](auto&& c_, auto&& predef_) {
                auto distr = ParameterDistributionConstant(c_);
                for (auto&& e : predef_) {
                    distr.add_predefined_sample(e);
                }
                return distr;
            },
            c, predef);
        if (p) {
            return success(p.value());
        }
        else {
            return p.as_failure();
        }
    }

    template <class IOContext>
    static IOResult<ParameterDistributionConstant> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("ParameterDistributionConstant");
        return deserialize_elements(io, obj);
    }

private:
    ScalarType m_constant;
};

template <class IOObj>
void details::SerializationVisitor<IOObj>::visit(const ParameterDistributionConstant& uniform_dist)
{
    obj.add_element("Type", std::string("Constant"));
    uniform_dist.serialize_elements(obj);
}

/**
 * deserialize a parameter distribution as a shared_ptr.
 * @see mio::deserialize
 */
template <class IOContext>
IOResult<std::shared_ptr<ParameterDistribution>> deserialize_internal(IOContext& io,
                                                                      Tag<std::shared_ptr<ParameterDistribution>>)
{
    auto obj  = io.expect_object("ParameterDistribution");
    auto type = obj.expect_element("Type", Tag<std::string>{});
    if (type) {
        if (type.value() == "Uniform") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionUniform::deserialize_elements(io, obj));
            return std::make_shared<ParameterDistributionUniform>(r);
        }
        else if (type.value() == "Normal") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionNormal::deserialize_elements(io, obj));
            return std::make_shared<ParameterDistributionNormal>(r);
        }
        else if (type.value() == "LogNormal") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionLogNormal::deserialize_elements(io, obj));
            return std::make_shared<ParameterDistributionLogNormal>(r);
        }
        else if (type.value() == "Exponential") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionExponential::deserialize_elements(io, obj));
            return std::make_shared<ParameterDistributionExponential>(r);
        }
        else if (type.value() == "Constant") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionConstant::deserialize_elements(io, obj));
            return std::make_shared<ParameterDistributionConstant>(r);
        }
        else {
            return failure(StatusCode::InvalidValue, "Type of ParameterDistribution " + type.value() + " not valid.");
        }
    }
    return failure(type.error());
}

} // namespace mio

#endif // PARAMETER_DISTRIBUTIONS_H
