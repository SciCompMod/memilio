#pragma once
#include <vector>
//#include <epidemiology/secir.h>
#include <epidemiology/eigen_util.h>

void save_result(std::vector<double> times, std::vector<Eigen::VectorXd> secir, std::string filename);

struct result
{
	std::vector<double> time;
	std::vector<std::vector<std::vector<double>>> groups;
	std::vector<std::vector<double>> total;
};

result read_result(std::string filename, int nb_groups);

