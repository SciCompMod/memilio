#pragma once
#include <vector>
//#include <epidemiology/secir.h>
#include <epidemiology/eigen_util.h>

void write_parameters(std::vector<epi::SecirParams> const& params, epi::ContactFrequencyMatrix const& cont_freq_matrix, double t0, double tmax, double dt, std::string filename);

struct file {
	int nb_groups;
	double t0;
	double tmax;
	double dt;
	std::vector<epi::SecirParams> const params;
	epi::ContactFrequencyMatrix const cont_freq_matrix;
};

file read_parameters(std::string filename);