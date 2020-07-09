#include <epidemiology/secir.h>
#include <epidemiology/logging.h>
#include <epidemiology/save_result.h>
#include <epidemiology/save_parameters.h>

int main()
{
    epi::set_log_level(epi::LogLevel::debug);
	
	file parameters = file{ read_parameters("Parameters.xml") };
	dist_params dists;

	

	double t0 = parameters.t0;
	double tmax = parameters.tmax;
	double dt = parameters.dt;
	int runs = parameters.runs;
	

	for (int run = 0; run < runs;++run)
	{
		std::vector<epi::SecirParams> params = parameters.params[run];
		epi::ContactFrequencyMatrix contact_freq_matrix = parameters.contact_freq_matrix[run];

		std::vector<Eigen::VectorXd> secir(0);
		std::vector<double> time = simulate(t0, tmax, dt, contact_freq_matrix, params, secir);


		save_result(time, secir, "Results/Result_run" + std::to_string(run) + ".h5");
		write_parameters(params, contact_freq_matrix, t0, tmax, dt, runs, "none", dists, "Results/Parameters_run" + std::to_string(run) + ".xml");
	}

}
