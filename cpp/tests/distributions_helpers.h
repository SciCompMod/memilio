#ifndef DISTRIBUTIONS_HELPERS_H
#define DISTRIBUTIONS_HELPERS_H

#include <epidemiology/parameter_studies/parameter_studies.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

void check_distribution(const epi::ParameterDistribution& dist, const epi::ParameterDistribution& dist_read);

//Mocks are not copyable because they need to store the call counters etc.
//ParameterDistribution must be copyable (for clone)
//so our mock consists of two classes
//the first class contains the actual mock logic but does not derive from anything
class MockParameterDistribution
{
public:
    MOCK_METHOD(double, get_rand_sample, (), ());
};
//the second class is clonable etc. and forwards calls to a stable instance of the first class
//templated so it allows StrictMock, NiceMock, etc.
//inherits from ParameterDistributionNormal instead of the base class so it is visitable
//all copies of an instance of this class will share the same mock and call counters etc.
template <class Mock = MockParameterDistribution>
class MockParameterDistributionRef : public epi::ParameterDistributionNormal
{
public:
    using epi::ParameterDistributionNormal::ParameterDistributionNormal;

    double get_rand_sample() override
    {
        return mock->get_rand_sample();
    }

    epi::ParameterDistribution* clone() const override
    {
        return new MockParameterDistributionRef(*this);
    }

    Mock& get_mock()
    {
        return *mock;
    }

private:
    std::shared_ptr<Mock> mock = std::make_shared<Mock>();
};

#endif