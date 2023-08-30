/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Martin J. Kuehn, Daniel Abele
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

#include "memilio/utils/logging.h"
#include "memilio/utils/visitor.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/io/io.h"

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
using ParameterDistributionVisitor = Visitor<class ParameterDistributionNormal, class ParameterDistributionUniform>;
using ConstParameterDistributionVisitor =
    ConstVisitor<class ParameterDistributionNormal, class ParameterDistributionUniform>;

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
    IOObj& obj;
};
} // namespace details

/*
 * Parameter Distribution class which contains the name of a variable as string
 * the lower bound and the upper bound as maximum admissible values and an enum
 * item with the name of the distribution
 */
class ParameterDistribution
{
public:
    ParameterDistribution(double lower_bound, double upper_bound)
        : m_lower_bound(lower_bound)
        , m_upper_bound(upper_bound)
    {
    }

    ParameterDistribution()
        : ParameterDistribution(0, 0)
    {
    }

    virtual ~ParameterDistribution() = default;

    void set_lower_bound(double lower_bound)
    {
        m_lower_bound = lower_bound;
    }

    void set_upper_bound(double upper_bound)
    {
        m_upper_bound = upper_bound;
    }

    void add_predefined_sample(double sample)
    {
        m_predefined_samples.push_back(sample);
    }

    void remove_predefined_samples()
    {
        m_predefined_samples.resize(0);
    }

    const std::vector<double>& get_predefined_samples() const
    {
        return m_predefined_samples;
    }

    double get_lower_bound() const
    {
        return m_lower_bound;
    }

    double get_upper_bound() const
    {
        return m_upper_bound;
    }

    /*
     * @brief returns a value for the given parameter distribution
     * in case some predefined samples are set, these values are taken
     * first, in case the vector of predefined values is empty, a 'real'
     * random sample is taken
     */
    double get_sample()
    {
        if (m_predefined_samples.size() > 0) {
            double rnumb = m_predefined_samples[0];
            m_predefined_samples.erase(m_predefined_samples.begin());
            return rnumb;
        }
        else {
            return get_rand_sample();
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

    virtual double get_rand_sample() = 0;

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
    double m_lower_bound; /*< A realistic lower bound on the given parameter */
    double m_upper_bound; /*< A realistic upper bound on the given parameter */
    std::vector<double>
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
    {
        m_mean         = 0;
        m_standard_dev = 1;
    }

    ParameterDistributionNormal(double mean, double standard_dev)
        : VisitableParameterDistribution<ParameterDistributionNormal>()
    {
        m_mean         = mean;
        m_standard_dev = standard_dev;
        check_quantiles(m_mean, m_standard_dev);
    }

    ParameterDistributionNormal(double lower_bound, double upper_bound, double mean)
        : VisitableParameterDistribution<ParameterDistributionNormal>(lower_bound, upper_bound)
    {
        m_mean         = mean;
        m_standard_dev = upper_bound; // set as to high and adapt then
        adapt_standard_dev(m_standard_dev);
    }

    ParameterDistributionNormal(double lower_bound, double upper_bound, double mean, double standard_dev)
        : VisitableParameterDistribution<ParameterDistributionNormal>(lower_bound, upper_bound)
    {
        m_mean         = mean;
        m_standard_dev = standard_dev;
        check_quantiles(m_mean, m_standard_dev);
    }

    void set_mean(double mean)
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
    bool check_quantiles(double& mean, double& standard_dev)
    {
        bool changed = false;

        changed = adapt_mean(mean);

        changed = adapt_standard_dev(standard_dev);

        if (changed && m_log_stddev_change) {
            log_warning("Standard deviation reduced to fit 99% of the distribution within [lowerbound,upperbound].");
        }

        return changed;
    }

    bool adapt_mean(double& mean)
    {
        bool changed = false;
        if (mean < m_lower_bound || mean > m_upper_bound) {
            mean    = 0.5 * (m_upper_bound - m_lower_bound);
            changed = true;
        }
        return changed;
    }

    // ensure that 0.99 % of the distribution are within lower bound and upper bound
    bool adapt_standard_dev(double& standard_dev)
    {
        bool changed = false;
        if (m_mean + standard_dev * m_quantile > m_upper_bound) {
            standard_dev = (m_upper_bound - m_mean) / m_quantile;
            changed      = true;
        }
        if (m_mean - standard_dev * m_quantile < m_lower_bound) {
            standard_dev = (m_mean - m_lower_bound) / m_quantile;
            changed      = true;
        }

        return changed;
    }

    void set_standard_dev(double standard_dev)
    {
        m_standard_dev = standard_dev;
    }

    double get_mean() const
    {
        return m_mean;
    }

    double get_standard_dev() const
    {
        return m_standard_dev;
    }

    void log_stddev_changes(bool log_stddev_change)
    {
        m_log_stddev_change = log_stddev_change;
    }

    /*
     * @brief gets a sample of a normally distributed variable
     * before sampling, it is verified that at least 99% of the
     * density function lie in the interval defined by the boundaries
     * otherwise the normal distribution is adapted
     */
    double get_rand_sample() override
    {
        //If ub = lb, sampling can only be succesful if mean = lb and dev = 0.
        //But this degenerate normal distribution is not allowed by the c++ standard.
        if (m_upper_bound == m_lower_bound) {
            return m_lower_bound;
        }

        if (check_quantiles(m_mean, m_standard_dev) || m_distribution.mean() != m_mean ||
            m_distribution.stddev() != m_standard_dev) {
            m_distribution = std::normal_distribution<double>{m_mean, m_standard_dev};
        }

        int i        = 0;
        int retries  = 10;
        double rnumb = m_distribution(thread_local_rng());
        while ((rnumb > m_upper_bound || rnumb < m_lower_bound) && i < retries) {
            rnumb = m_distribution(thread_local_rng());
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

    template<class IOObject>
    void serialize_elements(IOObject& obj) const
    {
        obj.add_element("Mean", m_mean);
        obj.add_element("StandardDev", m_standard_dev);
        obj.add_element("LowerBound", m_lower_bound);
        obj.add_element("UpperBound", m_upper_bound);
        obj.add_list("PredefinedSamples", m_predefined_samples.begin(), m_predefined_samples.end());
    }

    template<class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("ParameterDistributionNormal");
        serialize_elements(obj);
    }

    template<class IOContext, class IOObject>
    static IOResult<ParameterDistributionNormal> deserialize_elements(IOContext& io, IOObject& obj)
    {
        auto m      = obj.expect_element("Mean", Tag<double>{});
        auto s      = obj.expect_element("StandardDev", Tag<double>{});
        auto lb     = obj.expect_element("LowerBound", Tag<double>{});
        auto ub     = obj.expect_element("UpperBound", Tag<double>{});
        auto predef = obj.expect_list("PredefinedSamples", Tag<double>{});
        auto p      = apply(
            io,
            [](auto&& lb_, auto&& ub_, auto&& m_, auto&& s_, auto&& predef_) {
                auto distr = ParameterDistributionNormal(lb_, ub_, m_, s_);
                for (auto&& e : predef_) {
                    distr.add_predefined_sample(e);
                }
                return distr;
            },
            lb, ub, m, s, predef);
        if (p) {
            return success(p.value());
        }
        else {
            return p.as_failure();
        }
    }

    template<class IOContext>
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
    double m_mean; // the mean value of the normal distribution
    double m_standard_dev; // the standard deviation of the normal distribution
    constexpr static double m_quantile = 2.5758; // 0.995 quartile
    std::normal_distribution<double> m_distribution;
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
    ParameterDistributionUniform()
        : VisitableParameterDistribution<ParameterDistributionUniform>()
    {
    }

    ParameterDistributionUniform(double lower_bound, double upper_bound)
        : VisitableParameterDistribution<ParameterDistributionUniform>(lower_bound, upper_bound)
    {
    }

    /*
     * @brief gets a sample of a uniformly distributed variable
     */
    double get_rand_sample() override
    {
        if (m_distribution.max() != m_upper_bound || m_distribution.min() != m_lower_bound) {
            m_distribution = std::uniform_real_distribution<double>{m_lower_bound, m_upper_bound};
        }

        return m_distribution(thread_local_rng());
    }

    ParameterDistribution* clone() const override
    {
        return new ParameterDistributionUniform(*this);
    }

    template<class IOObject>
    void serialize_elements(IOObject& obj) const
    {
        obj.add_element("LowerBound", m_lower_bound);
        obj.add_element("UpperBound", m_upper_bound);
        obj.add_list("PredefinedSamples", m_predefined_samples.begin(), m_predefined_samples.end());
    }

    template<class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("ParameterDistributionUniform");
        serialize_elements(obj);
    }

    template<class IOContext, class IOObject>
    static IOResult<ParameterDistributionUniform> deserialize_elements(IOContext& io, IOObject& obj)
    {
        auto lb     = obj.expect_element("LowerBound", Tag<double>{});
        auto ub     = obj.expect_element("UpperBound", Tag<double>{});
        auto predef = obj.expect_list("PredefinedSamples", Tag<double>{});
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

    template<class IOContext>
    static IOResult<ParameterDistributionUniform> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("ParameterDistributionUniform");
        return deserialize_elements(io, obj);
    }

private:
    std::uniform_real_distribution<double> m_distribution;
};

template <class IOObj>
void details::SerializationVisitor<IOObj>::visit(const ParameterDistributionUniform& uniform_dist)
{
    obj.add_element("Type", std::string("Uniform"));
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
            BOOST_OUTCOME_TRY(r, ParameterDistributionUniform::deserialize_elements(io, obj));
            return std::make_shared<ParameterDistributionUniform>(r);
        }
        else if (type.value() == "Normal") {
            BOOST_OUTCOME_TRY(r, ParameterDistributionNormal::deserialize_elements(io, obj));
            return std::make_shared<ParameterDistributionNormal>(r);
        }
        else {
            return failure(StatusCode::InvalidValue, "Type of ParameterDistribution " + type.value() + " not valid.");
        }
    }
    return failure(type.error());
}

} // namespace mio

#endif // PARAMETER_DISTRIBUTIONS_H
