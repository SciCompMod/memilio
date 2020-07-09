#pragma once
#include <vector>
//#include <epidemiology/secir.h>
#include <epidemiology/eigen_util.h>

namespace epi
{

struct dist_params 
{
    std::vector<double> tinc;
    std::vector<double> tinfmild;
    std::vector<double> tserint;
    std::vector<double> thosp2home;
    std::vector<double> thome2hosp;
    std::vector<double> thosp2icu;
    std::vector<double> ticu2home;
    std::vector<double> tinfasy;
    std::vector<double> ticu2death;

    std::vector<double> inf_cont;
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<double> rho;
    std::vector<double> theta;
    std::vector<double> delta;
};

struct file 
{
    int nb_groups;
    int runs;
    double t0;
    double tmax;
    double dt;
    std::vector<std::vector<epi::SecirParams>> params;
    std::vector<epi::ContactFrequencyMatrix> contact_freq_matrix;
};

void write_parameters(std::vector<epi::SecirParams> const& params, epi::ContactFrequencyMatrix const& cont_freq_matrix,
                      double t0, double tmax, double dt, int runs, const std::string& dist, const dist_params& dists,
                      const std::string& filename);

double draw_number(double* values, char* dist);

file read_parameters(const std::string& filename);

} // namespace epi