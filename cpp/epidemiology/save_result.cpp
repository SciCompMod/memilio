#include <epidemiology/secir.h>
#include <epidemiology/damping.h>
#include <vector>

#include <iostream>
#include <string>
#include <H5Cpp.h>
using namespace H5;



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
	
	auto total = new double[n_data][n_compart];
	for (int i = 0;i < n_data;++i)
		for (int j = 0; j < n_compart;++j)
			total[i][j] = 0;

	for (int group = 0; group < nb_groups+1; ++group)
	{
		auto dset = new double[n_data][n_compart];
		if (group < nb_groups)
		{
			for (size_t irow = 0; irow < secir.size(); ++irow)
			{
				for (size_t compart = 0; compart < n_compart; ++compart)
				{
					dset[irow][compart] = secir[irow][compart + group * 8];
					total[irow][compart] += secir[irow][compart + group * 8];
				}
			}
		}
		// Create the data space for the dataset.
		hsize_t dims[2];               // dataset dimensions
		dims[0] = n_data;
		dims[1] = n_compart;
		DataSpace *dataspace = new DataSpace(RANK, dims);

		if (group == nb_groups)
			names.push_back("Total");
		else
			names.push_back("Group" + std::to_string(group + 1));

		// Create the dataset.      
		const H5std_string DATASET_NAME(names[group]);
		DataSet *dataset = new DataSet(file.createDataSet(DATASET_NAME,
			PredType::NATIVE_DOUBLE, *dataspace));

		if (group == nb_groups)
			dataset->write(total, PredType::NATIVE_DOUBLE);
		else
			dataset->write(dset, PredType::NATIVE_DOUBLE);

		delete dataset;
		delete dataspace;
	

	}
}

struct result 
{
	std::vector<double> time;
	std::vector<std::vector<std::vector<double>>> groups;
	std::vector<std::vector<double>> total;
};

result read_result(std::string filename, int nb_groups)
{

	
	const H5std_string FILE_NAME(filename);
	const int nb_compart = 8;

	std::vector<double> time_vec;

	

	H5File file(FILE_NAME, H5F_ACC_RDONLY);
	H5std_string DATASET_NAME_TIME("Time");
	DataSet dataset_time = file.openDataSet(DATASET_NAME_TIME);

	DataSpace filespace_time = dataset_time.getSpace();

	int rank_time = filespace_time.getSimpleExtentNdims();

	
	hsize_t dims_time[1];        // dataset dimensions
	rank_time = filespace_time.getSimpleExtentDims(dims_time);

	DataSpace mspace1(rank_time, dims_time);
	

	double* time = new double[dims_time[0]];
	
	dataset_time.read(time, PredType::NATIVE_DOUBLE, mspace1, filespace_time);
	
	for (int i = 0; i < dims_time[0];++i)
		time_vec.push_back(time[i]);


	std::vector<std::vector<std::vector<double>>> groups;
	
	
	for (int i = 0;i < nb_groups;++i)
	{
		H5std_string DATASET_NAME_GROUP("Group" + std::to_string(i+1));
		DataSet dataset_group = file.openDataSet(DATASET_NAME_GROUP);

		DataSpace filespace_group = dataset_group.getSpace();

		int rank_group = filespace_group.getSimpleExtentNdims();

		
		hsize_t dims_group[2];        // dataset dimensions
		rank_group = filespace_group.getSimpleExtentDims(dims_group);

		DataSpace mspace2(rank_group, dims_group);

		auto group = new double[dims_group[0]][nb_compart];
		
		dataset_group.read(group, PredType::NATIVE_DOUBLE, mspace2, filespace_group);
		
		std::vector<std::vector<double>> group_temp1;
		for (int l = 0; l < dims_group[0]; ++l)
		{
			std::vector<double> group_temp2;
			for (int k = 0; k < dims_group[1]; ++k)
			{
				group_temp2.push_back(group[l][k]);
			}
				
			group_temp1.push_back(group_temp2);
		}
		groups.push_back(group_temp1);

	}
	
	H5std_string DATASET_NAME_TOTAL("Total");
	DataSet dataset_total = file.openDataSet(DATASET_NAME_TOTAL);

	DataSpace filespace = dataset_total.getSpace();

	int rank_total = filespace.getSimpleExtentNdims();

	
	hsize_t dims_total[2];        // dataset dimensions
	rank_total = filespace.getSimpleExtentDims(dims_total);

	DataSpace mspace3(rank_total, dims_total);

	auto total = new double[dims_total[0]][nb_compart];
	dataset_total.read(total, PredType::NATIVE_DOUBLE, mspace3, filespace);

	std::vector<std::vector<double>> total_vec;
	for (int l = 0; l < dims_total[0]; ++l)
	{
		std::vector<double> total_temp;
		for (int k = 0; k < dims_total[1]; ++k)
			total_temp.push_back(total[l][k]);
		total_vec.push_back(total_temp);
	}

	return { time_vec, groups, total_vec };
}
