#include "reproduction_nb_secir.h"

double get_reproduction_number(Eigen::Index timept, mio::TimeSeries<double> y, mio::osecir::Parameters &params){
     //First get susceptibles in different groups at time t

        mio::ContactMatrixGroup &contact_matrix = params.get<mio::osecir::ContactPatterns>();

        Eigen::VectorXd susceptibles_at_t((size_t)params.get_num_groups());

        for (size_t agegrp = 0; agegrp < (size_t)params.get_num_groups(); agegrp++) {
            susceptibles_at_t[agegrp] =
                y.get_value(timept)[(Eigen::Index)mio::osecir::InfectionState::Susceptible +
                                          (int)mio::osecir::InfectionState::Count * agegrp];
        }

         //get population-sizes in the agegroups N_j

        Eigen::VectorXd population_agegroups((size_t)params.get_num_groups());
        for (size_t agegrp = 0; agegrp < (size_t)params.get_num_groups(); agegrp++) {
            for (size_t infectionstate = 0; infectionstate < (size_t)mio::osecir::InfectionState::Count;
                 infectionstate++)
                population_agegroups[agegrp] +=
                    y.get_value(timept)[infectionstate + (int)mio::osecir::InfectionState::Count * agegrp];
        }

         mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>::InternalArrayType beta =
            params.get<mio::osecir::RiskOfInfectionFromSymptomatic>().array();

        //Contact matrix \phi_{i,j}

        //T_{I_j} = time in compartment infected = TimeInfectedSymptoms ?
        //T_{C_j} = time in compartment Carrier = IncubationTime ?
        //p_i = TransmissionProbabilityonContact
        //beta = model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>().array() ?

        Eigen::MatrixXd eigenvalue_block((size_t)params.get_num_groups(), (size_t)params.get_num_groups());
        //Initialize eigenvalue_block
        for (int i = 0; i < eigenvalue_block.rows(); i++) {
            for (int j = 0; j < eigenvalue_block.cols(); j++) {
                eigenvalue_block(i, j) = susceptibles_at_t[i] *
                                         params.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)i] *
                                         contact_matrix.get_matrix_at(timept)(i, j) * 1 / (population_agegroups[i]) *
                                         (params.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)j] +
                                          beta[j] * params.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)j]);
            }
        }

        Eigen::ComplexEigenSolver<Eigen::MatrixXd> ces;

        ces.compute(eigenvalue_block);
        const Eigen::VectorXcd tempvector = ces.eigenvalues();

        Eigen::VectorXd tempvector1;
        tempvector1.resize(tempvector.size());
        //Do this later with some iterator
        for (int i = 0; i < tempvector.size(); i++) {
            tempvector1[i] = std::abs(tempvector[i]);
        }
        return tempvector1.maxCoeff();
    }

Eigen::VectorXd get_reproduction_numbers(mio::TimeSeries<double>y, mio::osecir::Parameters &params){
    Eigen::VectorXd temp(y.get_num_time_points());
    for(int i = 0; i < y.get_num_time_points(); i++){
        temp[i] = get_reproduction_number((Eigen::Index)i, y, params);
    }
    return temp;
}