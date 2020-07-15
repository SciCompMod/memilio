#include <epidemiology/secir.h>
#include <epidemiology/logging.h>
#include <epidemiology/save_result.h>
#include <epidemiology/save_parameters.h>

int main()
{
    // epi::set_log_level(epi::LogLevel::debug);
	
	// epi::file parameters = epi::file{ epi::read_parameters("Parameters.xml") };
	// epi::dist_params dists;

	

	// double t0 = parameters.t0;
	// double tmax = parameters.tmax;
	// double dt = parameters.dt;
	// int runs = parameters.runs;
	

	// for (int run = 0; run < runs;++run)
	// {
	// 	std::vector<epi::SecirParams> params = parameters.params[run];
	// 	epi::ContactFrequencyMatrix contact_freq_matrix = parameters.contact_freq_matrix[run];

	// 	std::vector<Eigen::VectorXd> secir(0);
	// 	std::vector<double> time = simulate(t0, tmax, dt, contact_freq_matrix, params, secir);


	// 	epi::save_result(time, secir, "Results/Result_run" + std::to_string(run) + ".h5");
	// 	epi::write_parameters(params, contact_freq_matrix, t0, tmax, dt, runs, dists, "Results/Parameters_run" + std::to_string(run) + ".xml");
	// }

}
