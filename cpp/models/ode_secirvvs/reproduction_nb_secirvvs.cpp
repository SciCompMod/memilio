#include "reproduction_nb_secirvvs.h"
#include "Eigen/src/Core/Matrix.h"
#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/time_series.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/parameters.h"

namespace mio{

    using Params = ParameterSet<IncubationPeriod, InfectedNoSymptomsToSymptoms, InfectedNoSymptomsToRecovered,
                 InfectedSymptomsToRecovered, InfectedSymptomsToSevere, SevereToCritical, SevereToRecovered,
                 CriticalToDead, CriticalToRecovered, RecoveredToSusceptible, ViralLoadDistributions,
                 InfectivityDistributions, DetectInfection, MaskProtection>;

    Eigen::VectorXd get_reproduction_number(Eigen::Index timept, mio::TimeSeries<ScalarType> result, mio::osecirvvs::Model model, Eigen::Ref<const Eigen::VectorXd> pop){
        //Calculate Delta_i
        Eigen::VectorXd Delta(model.parameters.get_num_groups());

        ContactMatrixGroup const& contact_matrix = model.parameters.get<mio::osecirvvs::ContactPatterns>();

        double season_val =
                    (1 + model.parameters.get<mio::osecirvvs::Seasonality>() *
                             sin(3.141592653589793 * (std::fmod((model.parameters.get<mio::osecirvvs::StartDay>() + timept), 365.0) / 182.5 + 0.5)));

        double cont_freq_eff =
                    season_val * contact_matrix.get_matrix_at(timept)(static_cast<Eigen::Index>((size_t)i),
                                                                 static_cast<Eigen::Index>((size_t)j));

        for(Eigen::Index i = 0; (mio::AgeGroup)i < model.parameters.get_num_groups(); i++){

                auto riskFromInfectedNoSymptoms = smoother_cosine(
                test_and_trace_required, params.get<TestAndTraceCapacity>(), params.get<TestAndTraceCapacity>() * 2,
                params.get<RelativeTransmissionNoSymptoms>()[i], 1.0);
                
                    double ext_inf_force_dummy = cont_freq_eff * divNj *
                                             model.parameters.template get<mio::osecirvvs::TransmissionProbabilityOnContact>()[(AgeGroup)i] *
                                             (riskFromInfectedNoSymptoms * (pop[INSNj] + pop[INSPIj] + pop[INSIIj]) +
                                              riskFromInfectedSymptomatic * (pop[ISyNj] + pop[ISyPIj] + pop[ISyIIj]));
            Delta[i] = 
        }
    }
}