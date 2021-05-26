#include "epidemiology/secir/parameter_space.h"
#include "epidemiology/utils/parameter_distributions.h"
#include "epidemiology/secir/secir.h"

namespace epi
{

void set_params_distributions_normal(
    SecirModel& model, double t0,
    double tmax, double dev_rel)
{
    auto set_distribution = [dev_rel](UncertainValue& v, double min_val = 0.001){
        v.set_distribution( ParameterDistributionNormal(std::max(min_val,
                                                       (1 - dev_rel * 2.6) * v),
                                                       (1 + dev_rel * 2.6) * v,
                                                       v,
                                                       dev_rel * v));
    };


    set_distribution(model.parameters.get<Seasonality>(), 0.0);
    set_distribution(model.parameters.get<ICUCapacity>());
    set_distribution(model.parameters.get<TestAndTraceCapacity>());

    // populations
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        for (auto j = Index<InfectionState>(0); j < Index<InfectionState>(InfectionState::Count); j++) {

            // don't touch S and D
            if ( j == InfectionState::Susceptible || j == InfectionState::Dead) {
                continue;
            }


            // variably sized groups
            set_distribution(model.populations[{i, j}]);
        }
    }

    // times
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {

        set_distribution(model.parameters.get<IncubationTime>()[i]);
        set_distribution(model.parameters.get<SerialInterval>()[i]);
        set_distribution(model.parameters.get<InfectiousTimeMild>()[i]);
        set_distribution(model.parameters.get<HospitalizedToHomeTime>()[i]);
        set_distribution(model.parameters.get<HomeToHospitalizedTime>()[i]);
        set_distribution(model.parameters.get<InfectiousTimeAsymptomatic>()[i]);
        set_distribution(model.parameters.get<HospitalizedToICUTime>()[i]);
        set_distribution(model.parameters.get<ICUToHomeTime>()[i]);
        set_distribution(model.parameters.get<ICUToDeathTime>()[i]);

        set_distribution(model.parameters.get<InfectionProbabilityFromContact>()[i]);
        set_distribution(model.parameters.get<RelativeCarrierInfectability>()[i]);
        set_distribution(model.parameters.get<AsymptoticCasesPerInfectious>()[i]);
        set_distribution(model.parameters.get<RiskOfInfectionFromSympomatic>()[i]);
        set_distribution(model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[i]);
        set_distribution(model.parameters.get<DeathsPerHospitalized>()[i]);
        set_distribution(model.parameters.get<HospitalizedCasesPerInfectious>()[i]);
        set_distribution(model.parameters.get<ICUCasesPerHospitalized>()[i]);
    }

    // maximum number of dampings; to avoid overfitting only allow one damping for every 10 days simulated
    // damping base values are between 0.1 and 1; diagonal values vary lie in the range of 0.6 to 1.4 times the base value
    // off diagonal values vary between 0.7 to 1.1 of the corresponding diagonal value (symmetrization is conducted)
    model.parameters.get<ContactPatterns>().set_distribution_damp_nb(ParameterDistributionUniform(1, (tmax - t0) / 10));
    model.parameters.get<ContactPatterns>().set_distribution_damp_days(ParameterDistributionUniform(t0, tmax));
    model.parameters.get<ContactPatterns>().set_distribution_damp_diag_base(ParameterDistributionUniform(0.0, 0.9));
    model.parameters.get<ContactPatterns>().set_distribution_damp_diag_rel(ParameterDistributionUniform(0.0, 0.4));
    model.parameters.get<ContactPatterns>().set_distribution_damp_offdiag_rel(ParameterDistributionUniform(0.0, 0.3));
}


void draw_sample_demographics(SecirModel& model)
{
    model.parameters.get<ICUCapacity>().draw_sample();
    model.parameters.get<TestAndTraceCapacity>().draw_sample();

    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        double group_total = model.populations.get_group_total(i);

        model.populations[{i, InfectionState::Exposed}].draw_sample();
        model.populations[{i, InfectionState::Carrier}].draw_sample();
        model.populations[{i, InfectionState::Infected}].draw_sample();
        model.populations[{i, InfectionState::Hospitalized}].draw_sample();
        model.populations[{i, InfectionState::ICU}].draw_sample();
        model.populations[{i, InfectionState::Recovered}].draw_sample();

        // no sampling for dead and total numbers
        // [...]

        model.populations.set_difference_from_group_total<AgeGroup>({i, InfectionState::Susceptible}, group_total);
        model.populations.set_difference_from_group_total<AgeGroup>({i, InfectionState::Susceptible},
                                                                    model.populations.get_group_total(i));
    }
}

void draw_sample_infection(SecirModel& model)
{
    model.parameters.get<Seasonality>().draw_sample();

    //not age dependent
    model.parameters.get<IncubationTime>()[AgeGroup(0)].draw_sample();
    model.parameters.get<SerialInterval>()[AgeGroup(0)].draw_sample();
    model.parameters.get<InfectiousTimeMild>()[AgeGroup(0)].draw_sample();
    model.parameters.get<HospitalizedToICUTime>()[AgeGroup(0)].draw_sample();
    model.parameters.get<RelativeCarrierInfectability>()[AgeGroup(0)].draw_sample();
    model.parameters.get<RiskOfInfectionFromSympomatic>()[AgeGroup(0)].draw_sample();
    model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[AgeGroup(0)].draw_sample();

    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        //not age dependent
        model.parameters.get<IncubationTime>()[i] =
            model.parameters.get<IncubationTime>()[AgeGroup(0)];
        model.parameters.get<SerialInterval>()[i] =
            model.parameters.get<SerialInterval>()[AgeGroup(0)];
        model.parameters.get<InfectiousTimeMild>()[i] =
            model.parameters.get<InfectiousTimeMild>()[AgeGroup(0)];
        model.parameters.get<HospitalizedToICUTime>()[i] =
            model.parameters.get<HospitalizedToICUTime>()[AgeGroup(0)];
        model.parameters.get<RelativeCarrierInfectability>()[i] =
            model.parameters.get<RelativeCarrierInfectability>()[AgeGroup(0)];
        model.parameters.get<RiskOfInfectionFromSympomatic>()[i] =
            model.parameters.get<RiskOfInfectionFromSympomatic>()[AgeGroup(0)];
        model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[i] =
            model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[AgeGroup(0)];

        //age dependent
        model.parameters.get<HospitalizedToHomeTime>()[i].draw_sample(); // here: home=recovered
        model.parameters.get<HomeToHospitalizedTime>()[i].draw_sample(); // here: home=infectious
        model.parameters.get<InfectiousTimeAsymptomatic>()[i].draw_sample();
        model.parameters.get<ICUToDeathTime>()[i].draw_sample();
        model.parameters.get<ICUToHomeTime>()[i].draw_sample();

        model.parameters.get<InfectionProbabilityFromContact>()[i].draw_sample();
        model.parameters.get<AsymptoticCasesPerInfectious>()[i].draw_sample();
        model.parameters.get<DeathsPerHospitalized>()[i].draw_sample();
        model.parameters.get<HospitalizedCasesPerInfectious>()[i].draw_sample();
        model.parameters.get<ICUCasesPerHospitalized>()[i].draw_sample();
    }
}

void draw_sample(SecirModel& model)
{
    draw_sample_infection(model);
    draw_sample_demographics(model);
    model.parameters.get<ContactPatterns>().draw_sample();
    model.apply_constraints();
}

} // namespace epi
