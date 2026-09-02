#include "parameter_setter.h"
#include "constants.h"

#include "memilio/utils/parameter_distributions.h"

#include <cassert>
#include <cmath>

std::pair<double, double> get_mu_and_sigma(std::pair<double, double> mean_and_std)
{
    const auto mean    = mean_and_std.first;
    const auto stddev  = mean_and_std.second;
    const double mu    = std::log(mean * mean / std::sqrt(mean * mean + stddev * stddev));
    const double sigma = std::sqrt(std::log(1 + stddev * stddev / (mean * mean)));
    return {mu, sigma};
}

namespace
{
/// @brief Assign a log normal distribution given by the mean and standard deviation of the values.
template <class Parameter>
void set_lognormal(mio::abm::Parameters& params, double mean, double stddev)
{
    const auto mu_and_sigma = get_mu_and_sigma({mean, stddev});
    params.get<Parameter>() = mio::ParameterDistributionLogNormal(mu_and_sigma.first, mu_and_sigma.second);
}
} // namespace

void set_parameters(mio::abm::Parameters& params)
{
    // Duration of the infection states. What the panvXabmSim branch called IncubationPeriod is
    // called TimeExposedToNoSymptoms on main.
    set_lognormal<mio::abm::TimeExposedToNoSymptoms>(params, 4.5, 1.5);
    set_lognormal<mio::abm::TimeInfectedNoSymptomsToSymptoms>(params, 1.1, 0.9);
    set_lognormal<mio::abm::TimeInfectedNoSymptomsToRecovered>(params, 8.0, 2.0);
    set_lognormal<mio::abm::TimeInfectedSymptomsToSevere>(params, 6.6, 4.9);
    set_lognormal<mio::abm::TimeInfectedSymptomsToRecovered>(params, 8.0, 2.0);
    set_lognormal<mio::abm::TimeInfectedSevereToCritical>(params, 1.5, 2.0);
    set_lognormal<mio::abm::TimeInfectedSevereToRecovered>(params, 18.1, 6.3);
    set_lognormal<mio::abm::TimeInfectedCriticalToDead>(params, 10.7, 4.8);
    set_lognormal<mio::abm::TimeInfectedCriticalToRecovered>(params, 18.1, 6.3);

    // Transition probabilities per age group.
    constexpr auto virus = mio::abm::VirusVariant::Wildtype;

    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{virus, age_group_0_to_4}]   = 0.50;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{virus, age_group_5_to_14}]  = 0.55;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{virus, age_group_15_to_34}] = 0.60;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{virus, age_group_35_to_59}] = 0.70;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{virus, age_group_60_to_79}] = 0.83;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{virus, age_group_80_plus}]  = 0.90;

    params.get<mio::abm::SeverePerInfectedSymptoms>()[{virus, age_group_0_to_4}]   = 0.02;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{virus, age_group_5_to_14}]  = 0.03;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{virus, age_group_15_to_34}] = 0.04;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{virus, age_group_35_to_59}] = 0.07;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{virus, age_group_60_to_79}] = 0.17;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{virus, age_group_80_plus}]  = 0.24;

    params.get<mio::abm::CriticalPerInfectedSevere>()[{virus, age_group_0_to_4}]   = 0.10;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{virus, age_group_5_to_14}]  = 0.11;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{virus, age_group_15_to_34}] = 0.12;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{virus, age_group_35_to_59}] = 0.14;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{virus, age_group_60_to_79}] = 0.33;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{virus, age_group_80_plus}]  = 0.62;

    params.get<mio::abm::DeathsPerInfectedCritical>()[{virus, age_group_0_to_4}]   = 0.12;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{virus, age_group_5_to_14}]  = 0.13;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{virus, age_group_15_to_34}] = 0.15;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{virus, age_group_35_to_59}] = 0.26;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{virus, age_group_60_to_79}] = 0.40;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{virus, age_group_80_plus}]  = 0.48;

    // Transmission. The scale of the viral shed is overwritten per run from the --infection_k option.
    params.get<mio::abm::InfectionRateFromViralShed>()[{virus}] = 1.0;
    params.get<mio::abm::AerosolTransmissionRates>()            = 0.0;

    // Mobility. These have to be set explicitly, the defaults on main let nobody go to school or work.
    // The 15-34 group appears in both sets because CityBuilder assigns it to schools when there are
    // not enough 5-14 year olds to fill them; each Person is assigned either a School or a Work.
    params.get<mio::abm::AgeGroupGotoSchool>()                     = false;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14]  = true;
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_15_to_34] = true;
    params.get<mio::abm::AgeGroupGotoWork>()                       = false;
    params.get<mio::abm::AgeGroupGotoWork>().set_multiple(
        {age_group_15_to_34, age_group_35_to_59, age_group_60_to_79, age_group_80_plus}, true);
}

namespace
{

/// @brief The four baseline contact matrices a region is described by.
struct ContactMatrixSet {
    Eigen::MatrixXd home;
    Eigen::MatrixXd school;
    Eigen::MatrixXd work;
    Eigen::MatrixXd other;
};

/// @brief Baseline contact matrices for Germany, indexed (receiving age group, transmitting age group).
ContactMatrixSet contacts_germany(Eigen::Index n)
{
    ContactMatrixSet c{Eigen::MatrixXd::Zero(n, n), Eigen::MatrixXd::Zero(n, n),
                       Eigen::MatrixXd::Zero(n, n), Eigen::MatrixXd::Zero(n, n)};

    c.home.row(0) << 0.4413, 0.0504, 1.2383, 0.8033, 0.0494, 0.0017;
    c.home.row(1) << 0.0485, 0.7616, 0.6532, 1.1614, 0.0256, 0.0013;
    c.home.row(2) << 0.1800, 0.1795, 0.8806, 0.6413, 0.0429, 0.0032;
    c.home.row(3) << 0.0495, 0.2639, 0.5189, 0.8277, 0.0679, 0.0014;
    c.home.row(4) << 0.0087, 0.0394, 0.1417, 0.3834, 0.7064, 0.0447;
    c.home.row(5) << 0.0292, 0.0648, 0.1248, 0.4179, 0.3497, 0.1544;

    c.school.row(0) << 1.1165, 0.2741, 0.2235, 0.1028, 0.0007, 0.0000;
    c.school.row(1) << 0.1627, 1.9412, 0.2431, 0.1780, 0.0130, 0.0000;
    c.school.row(2) << 0.0148, 0.1646, 1.1266, 0.0923, 0.0074, 0.0000;
    c.school.row(3) << 0.0367, 0.1843, 0.3265, 0.0502, 0.0021, 0.0005;
    c.school.row(4) << 0.0004, 0.0370, 0.0115, 0.0014, 0.0039, 0.0000;
    c.school.row(5) << 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000;

    c.work.row(0) << 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000;
    c.work.row(1) << 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000;
    c.work.row(2) << 0.0000, 0.0127, 1.7570, 1.6050, 0.0133, 0.0000;
    c.work.row(3) << 0.0000, 0.0020, 1.0311, 2.3166, 0.0098, 0.0000;
    c.work.row(4) << 0.0000, 0.0002, 0.0194, 0.0325, 0.0003, 0.0000;
    c.work.row(5) << 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000;

    c.other.row(0) << 0.5170, 0.3997, 0.7957, 0.9958, 0.3239, 0.0428;
    c.other.row(1) << 0.0632, 0.9121, 0.3254, 0.4731, 0.2355, 0.0148;
    c.other.row(2) << 0.0336, 0.1604, 1.7529, 0.8622, 0.1440, 0.0077;
    c.other.row(3) << 0.0204, 0.1444, 0.5738, 1.2127, 0.3433, 0.0178;
    c.other.row(4) << 0.0371, 0.0393, 0.4171, 0.9666, 0.7495, 0.0257;
    c.other.row(5) << 0.0791, 0.0800, 0.3480, 0.5588, 0.2769, 0.0180;

    return c;
}

/// @brief Baseline contact matrices for France, indexed (receiving age group, transmitting age group).
ContactMatrixSet contacts_france(Eigen::Index n)
{
    ContactMatrixSet c{Eigen::MatrixXd::Zero(n, n), Eigen::MatrixXd::Zero(n, n),
                       Eigen::MatrixXd::Zero(n, n), Eigen::MatrixXd::Zero(n, n)};

    c.home.row(0) << 0.6881, 0.6771, 1.2965, 0.9261, 0.0337, 0.0034;
    c.home.row(1) << 0.2257, 1.6804, 0.7570, 1.4088, 0.0235, 0.0022;
    c.home.row(2) << 0.2563, 0.3517, 1.4941, 0.7716, 0.0381, 0.0020;
    c.home.row(3) << 0.2096, 0.6996, 0.8293, 1.2112, 0.0640, 0.0075;
    c.home.row(4) << 0.1710, 0.4608, 0.4983, 0.6181, 0.8598, 0.0665;
    c.home.row(5) << 0.1990, 0.6331, 0.5450, 1.1572, 0.4189, 0.2800;

    c.school.row(0) << 2.4233, 0.3734, 0.3248, 0.3342, 0.0032, 0.0000;
    c.school.row(1) << 0.1972, 3.3900, 0.1830, 0.2888, 0.0054, 0.0001;
    c.school.row(2) << 0.0317, 0.4252, 1.3594, 0.1850, 0.0051, 0.0002;
    c.school.row(3) << 0.1477, 0.4910, 0.4562, 0.2218, 0.0061, 0.0001;
    c.school.row(4) << 0.0207, 0.0356, 0.0564, 0.0603, 0.0277, 0.0035;
    c.school.row(5) << 0.0000, 0.0212, 0.0289, 0.0000, 0.0000, 0.0000;

    c.work.row(0) << 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000;
    c.work.row(1) << 0.0000, 0.0074, 0.0222, 0.0254, 0.0000, 0.0000;
    c.work.row(2) << 0.0000, 0.0235, 2.7163, 2.4826, 0.0050, 0.0000;
    c.work.row(3) << 0.0000, 0.0249, 2.0297, 3.5135, 0.0052, 0.0000;
    c.work.row(4) << 0.0000, 0.0005, 0.0133, 0.0262, 0.0001, 0.0000;
    c.work.row(5) << 0.0000, 0.0000, 0.0001, 0.0001, 0.0000, 0.0000;

    c.other.row(0) << 0.6935, 0.4678, 0.9381, 1.1670, 0.3841, 0.0317;
    c.other.row(1) << 0.2217, 2.2229, 0.8144, 1.0773, 0.2950, 0.0378;
    c.other.row(2) << 0.1051, 0.4410, 2.7757, 1.2325, 0.1793, 0.0230;
    c.other.row(3) << 0.0847, 0.2358, 1.1066, 1.8774, 0.4578, 0.0398;
    c.other.row(4) << 0.0495, 0.1284, 0.7258, 1.5573, 1.0805, 0.1023;
    c.other.row(5) << 0.0446, 0.1150, 0.4335, 1.0023, 0.8198, 0.1387;

    return c;
}

/// @brief Baseline contact matrices for the United States, indexed (receiving age group, transmitting age group).
ContactMatrixSet contacts_united_states(Eigen::Index n)
{
    ContactMatrixSet c{Eigen::MatrixXd::Zero(n, n), Eigen::MatrixXd::Zero(n, n),
                       Eigen::MatrixXd::Zero(n, n), Eigen::MatrixXd::Zero(n, n)};

    c.home.row(0) << 0.6197, 0.7925, 1.1687, 0.9485, 0.0335, 0.0040;
    c.home.row(1) << 0.2460, 1.7677, 0.8041, 1.3347, 0.0308, 0.0032;
    c.home.row(2) << 0.2832, 0.4564, 1.5297, 0.7018, 0.0449, 0.0028;
    c.home.row(3) << 0.2613, 0.8550, 0.9010, 1.1937, 0.0834, 0.0095;
    c.home.row(4) << 0.2117, 0.5553, 0.5944, 0.7254, 0.6861, 0.0641;
    c.home.row(5) << 0.1787, 0.5972, 0.5194, 0.9997, 0.2481, 0.2714;

    c.school.row(0) << 1.1966, 0.2696, 0.2404, 0.2202, 0.0036, 0.0000;
    c.school.row(1) << 0.1381, 3.9843, 0.2120, 0.2775, 0.0077, 0.0001;
    c.school.row(2) << 0.0240, 0.5101, 1.8745, 0.1804, 0.0066, 0.0003;
    c.school.row(3) << 0.0982, 0.4768, 0.4823, 0.1892, 0.0071, 0.0001;
    c.school.row(4) << 0.0251, 0.0527, 0.0906, 0.0676, 0.0291, 0.0044;
    c.school.row(5) << 0.0000, 0.0211, 0.0289, 0.0000, 0.0000, 0.0000;

    c.work.row(0) << 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000;
    c.work.row(1) << 0.0000, 0.0402, 0.0502, 0.0548, 0.0000, 0.0000;
    c.work.row(2) << 0.0000, 0.0526, 2.5209, 2.3056, 0.0387, 0.0000;
    c.work.row(3) << 0.0000, 0.0664, 1.9258, 3.4223, 0.0442, 0.0000;
    c.work.row(4) << 0.0000, 0.0105, 0.1102, 0.2350, 0.0047, 0.0000;
    c.work.row(5) << 0.0000, 0.0000, 0.0001, 0.0001, 0.0000, 0.0000;

    c.other.row(0) << 0.7819, 0.5386, 1.0918, 1.2430, 0.3503, 0.0309;
    c.other.row(1) << 0.2522, 2.6569, 0.9886, 1.1725, 0.2698, 0.0380;
    c.other.row(2) << 0.1245, 0.5258, 3.4633, 1.3746, 0.1660, 0.0232;
    c.other.row(3) << 0.0913, 0.2606, 1.2291, 1.9120, 0.3962, 0.0368;
    c.other.row(4) << 0.0469, 0.1202, 0.7026, 1.3731, 0.7760, 0.0768;
    c.other.row(5) << 0.0315, 0.0835, 0.3139, 0.6754, 0.4473, 0.0847;

    return c;
}

ContactMatrixSet contacts_for(CityParameters::Region region, Eigen::Index n)
{
    switch (region) {
    case CityParameters::Region::France:
        return contacts_france(n);
    case CityParameters::Region::UnitedStates:
        return contacts_united_states(n);
    case CityParameters::Region::Germany:
        break;
    }
    return contacts_germany(n);
}

} // namespace

void set_local_parameters(mio::abm::Model& model, CityParameters::Region region)
{
    const auto n = static_cast<Eigen::Index>((size_t)model.parameters.get_num_groups());
    assert(n == static_cast<Eigen::Index>(num_age_groups) &&
           "The contact matrices of this simulation are defined for six age groups.");

    const auto contacts = contacts_for(region, n);
    // Location types that are not listed below fall back to one contact per age group pair and day.
    const Eigen::MatrixXd contacts_default = Eigen::MatrixXd::Constant(n, n, 1.0);

    for (auto& loc : model.get_locations()) {
        auto& baseline = loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline();
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            // 17 hours per day, scaled by the contact intensity at home.
            baseline = contacts.home * 1.4 * 15.0;
            break;
        case mio::abm::LocationType::School:
            // 5 hours per day.
            baseline = contacts.school * 4.8;
            break;
        case mio::abm::LocationType::Work:
            // 7 hours per day.
            baseline = contacts.work * 3.0 * 0.5;
            break;
        case mio::abm::LocationType::SocialEvent:
            // 4 hours per day, scaled by the share of "other" contacts and the contact intensity.
            baseline = contacts.other * 1.2 * 2.0 * 6.0;
            break;
        case mio::abm::LocationType::BasicsShop:
            // 2 hours per day, scaled by the share of "other" contacts and the contact intensity.
            baseline = contacts.other * 0.8 * 0.33 * 12.0;
            break;
        default:
            baseline = contacts_default;
            break;
        }
    }
}
