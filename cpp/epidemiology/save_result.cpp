#include <epidemiology/secir.h>
#include <epidemiology/damping.h>
#include <vector>

#include <iostream>
#include <string>
#include <H5Cpp.h>
using namespace H5;

#include <tixi.h>


void save_result(std::vector<double> times, std::vector<Eigen::VectorXd> secir, std::string filename)
{
	const H5std_string	FILE_NAME(filename);
	const int	 RANK = 2;

	const int n_data = secir.size();
	const int n_compart = 8;
	const int nb_groups = secir[0].size() / n_compart;

	std::vector<std::string> names;

	// Create a new file using the default property lists. 
	H5File file(FILE_NAME, H5F_ACC_TRUNC);

	hsize_t dims[1];               // dataset dimensions
	dims[0] = n_data;
	DataSpace *dataspace = new DataSpace(1, dims);

	const H5std_string DATASET_NAME("Time");
	DataSet *dataset = new DataSet(file.createDataSet(DATASET_NAME,
		PredType::NATIVE_DOUBLE, *dataspace));

	auto time = new double[n_data];
	for (size_t irow = 0; irow < secir.size(); ++irow)
	{
		time[irow] = times[irow];
	}

	dataset->write(time, PredType::NATIVE_DOUBLE);
	

	for (int group = 0; group < nb_groups; ++group)
	{
		auto dset = new double[n_data][n_compart];
		names.push_back("Group" + std::to_string(group + 1));
		for (size_t irow = 0; irow < secir.size(); ++irow)
		{
			dset[irow][0] = secir[irow][0 + group * 8];
			dset[irow][1] = secir[irow][1 + group * 8];
			dset[irow][2] = secir[irow][2 + group * 8];
			dset[irow][3] = secir[irow][3 + group * 8];
			dset[irow][4] = secir[irow][4 + group * 8];
			dset[irow][5] = secir[irow][5 + group * 8];
			dset[irow][6] = secir[irow][6 + group * 8];
			dset[irow][7] = secir[irow][7 + group * 8];
		}

		// Create the data space for the dataset.
		hsize_t dims[2];               // dataset dimensions
		dims[0] = n_data;
		dims[1] = n_compart;
		DataSpace *dataspace = new DataSpace(RANK, dims);

		// Create the dataset.      
		const H5std_string DATASET_NAME(names[group]);
		DataSet *dataset = new DataSet(file.createDataSet(DATASET_NAME,
			PredType::NATIVE_DOUBLE, *dataspace));

		dataset->write(dset, PredType::NATIVE_DOUBLE);


		delete dataset;
		delete dataspace;

	}
}

