#pragma once
#include <vector>
//#include <epidemiology/secir.h>
#include <epidemiology/eigen_util.h>

void save_result(std::vector<double> times, std::vector<Eigen::VectorXd> secir, std::string filename);

