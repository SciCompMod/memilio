#pragma once
#include <vector>
//#include <epidemiology/secir.h>
#include <epidemiology/eigen_util.h>

namespace epi
{

struct dist_params 
{
	std::string dist_total = "none";
	std::string dist_exposed = "none";
	std::string dist_carrier = "none";
	std::string dist_infectious = "none";
	std::string dist_hospital = "none";
	std::string dist_icu = "none";
	std::string dist_recovered = "none";
	std::string dist_dead = "none";
	
	std::string dist_tinc = "none";
	std::string dist_tinfmild = "none";
	std::string dist_tserint = "none";
	std::string dist_thosp2home = "none";
	std::string dist_thome2hosp = "none";
	std::string dist_thosp2icu = "none";
	std::string dist_ticu2home = "none";
	std::string dist_tinfasy = "none";
	std::string dist_ticu2death = "none";
	
	std::string dist_inf_cont = "none";
	std::string dist_alpha = "none";
	std::string dist_beta = "none";
	std::string dist_rho = "none";
	std::string dist_theta = "none";
	std::string dist_delta = "none";
	
	std::vector<double> total;
	std::vector<double> exposed;
	std::vector<double> carrier;
	std::vector<double> infectious;
	std::vector<double> hospital;
	std::vector<double> icu;
	std::vector<double> recovered;
	std::vector<double> dead;
	
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

ContactFrequencyVariableElement read_contact(const std::string& filename);
ParameterDistribution read_dist(const std::string& filename, std::string& path);
file read_global_params(const std::string& filename);

void write_parameters(std::vector<epi::SecirParams> const& params, epi::ContactFrequencyMatrix const& cont_freq_matrix,
                      double t0, double tmax, double dt, int runs, const std::string& dist, const dist_params& dists,
                      const std::string& filename);

double draw_number(double* values, char* dist);

file read_parameters(const std::string& filename);

} // namespace epi