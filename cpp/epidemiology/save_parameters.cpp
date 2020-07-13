#include <epidemiology/secir.h>
#include <epidemiology/damping.h>
#include <epidemiology/parameter_studies/parameter_space.h>
#include <vector>

#include <iostream>
#include <string>
#include <random>

#include <tixi.h>

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

void write_dist(const TixiDocumentHandle& handle, const std::string& path, const std::string& name, const std::string& dist, double param, const std::vector<double>& dist_param){

	tixiCreateElement(handle, path.c_str(), name.c_str());
	std::string full_path = path + "/" + name
	tixiAddTextElement(handle, (full_path).c_str(), "Distribution", dist.c_str());
	if (strcmp("none", dist) == 0){
		tixiAddDoubleElement(handle, full_path.c_str(), "Mean", param, "%g")
	}
	else if (strcmp("normal", dist) == 0){
		tixiAddDoubleElement(handle, full_path.c_str(), "Mean", param, "%g");
		tixiAddDoubleElement(handle, full_path.c_str(), "Deviation", dist_param[2], "%g");
		tixiAddDoubleElement(handle, full_path.c_str(), "Min", dist_param[0], "%g");
		tixiAddDoubleElement(handle, full_path.c_str(), "Max", dist_param[1], "%g");
	}
	else if (strcmp("uniform", dist) == 0){
		tixiAddDoubleElement(handle, full_path.c_str(), "Min", dist_param[0], "%g");
		tixiAddDoubleElement(handle, full_path.c_str(), "Max", dist_param[1], "%g");
	}
}	

struct file
{
	int nb_groups;
	int runs;
	double t0;
	double tmax;
	double dt;
};


file read_global_params(const std::string& filename){
	TixiDocumentHandle handle = -1;
	tixiOpenDocument(filename.c_str(), &handle);
	int nb_groups;
	int runs;
	double t0;
	double tmax;
	double dt;


	tixiGetIntegerElement(handle, "/Parameters/NumberOfGroups", &nb_groups);
	tixiGetIntegerElement(handle, "/Parameters/Runs", &runs);
	tixiGetTextElement(handle, "/Parameters/Distribution", &dist);
	tixiGetDoubleElement(handle, "/Parameters/T0", &t0);
	tixiGetDoubleElement(handle, "/Parameters/TMax", &tmax);
	tixiGetDoubleElement(handle, "/Parameters/dt", &dt);
	return {nb_groups, runs, t0, tmax, dt};
}


ParameterDistribution read_dist(const std::string& filename, std::string& path){
	
	TixiDocumentHandle handle = -1;
	tixiOpenDocument(filename.c_str(), &handle);
	char* dist;
	tixiGetTextElement(handle, (path + "/Distribution").c_str(), &dist);
	
	if (strcmp("none", dist) == 0){
		double mean;
		tixiGetDoubleElement(handle, (path + "/Mean").c_str(), &mean);
		Distribution = std::make_unique<ParameterDistributionNormal>(
            ParameterDistributionNormal(mean, 0));
	}
	else if (strcmp("normal", dist) == 0){
		double mean;
		double dev;
		double min;
		double max;
		tixiGetDoubleElement(handle, (path + "/Mean").c_str(), &mean);
		tixiGetDoubleElement(handle, (path + "/Deviation").c_str(), &dev);
		tixiGetDoubleElement(handle, (path + "/Min").c_str(), &min);
		tixiGetDoubleElement(handle, (path + "/Max").c_str(), &max);
		Distribution = std::make_unique<ParameterDistributionNormal>(
            ParameterDistributionNormal(min, max, mean, dev));
	}
	else if (strcmp("uniform", dist) == 0){
		double min;
		double max;
		tixiGetDoubleElement(handle, (path + "/Min").c_str(), &min);
		tixiGetDoubleElement(handle, (path + "/Max").c_str(), &max);
		Distribution = std::make_unique<ParameterDistributionUniform>(
            ParameterDistributionUniform(min, max));
	}
	return Distribution;
	
}

ContactFrequencyVariableElement read_contact(const std::string& filename){
	epi::ContactFrequencyMatrix contact_freq_matrix{ (size_t)nb_groups };
	
	for (size_t i = 0; i < nb_groups; i++) {
		double* contact = new double[nb_groups];
		tixiGetFloatVector(handle, ("/Parameters/ContactRateGroup" + std::to_string(i + 1)).c_str(), &contact, nb_groups);
		for (int j = 0; j < nb_groups; ++j)
		{
			contact_freq_matrix.set_cont_freq(contact[j], i, j);
		}
		delete contact;
	}
	
	
	TixiDocumentHandle handle = -1;
	tixiOpenDocument(filename.c_str(), &handle);
	int nb_damp;
	tixiGetIntegerElement(handle, "/Parameters/NumberOfDampings", &nb_damp);
	std::unique_ptr<ContactFrequencyVariableElement> m_cont_freq_matrix_variable;
	
	nb_damp_dist = std::make_unique<ParameterDistributionUniform>(ParameterDistributionUniform(1, (tmax - t0) / 10));
	 day_dist = std::make_unique<ParameterDistributionUniform>(ParameterDistributionUniform(t0, tmax));
	 damp_diag_base_dist = std::make_unique<ParameterDistributionUniform>(ParameterDistributionUniform(0.1, 1));
	 damp_diag_rel_dist = std::make_unique<ParameterDistributionUniform>(ParameterDistributionUniform(0.6, 1.4));
	 damp_offdiag_rel_dist = std::make_unique<ParameterDistributionUniform>(ParameterDistributionUniform(0.7, 1.1));
	 
	 nb_damp_dist.add_predefined_sample(nb_damp);
	 for (int i = 0; i < nb_damp; ++i)
	 {
		double day;
		double damping;

		tixiGetDoubleElement(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1) + "/Day").c_str(), &day);
		day_dist.add_predefined_sample(day);
		for (int k = 0; k < nb_groups; ++k)
		{
			tixiGetDoubleElement(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1) + "/Group" + std::to_string(k + 1)).c_str(), &damping);
			damp_diag_base.add_predefined_sample(damping);
			damp_diag_rel_dist.add_predefined_sample(1.0);
			for (int l = k+1; l<nb_groups; ++l)
				damp_offdiag_rel_dist.add_predefined_sample(1.0);
		}
	 }
	 
	 
	// maximum number of dampings; to avoid overfitting only allow one damping for every 10 days simulated
    // damping base values are between 0.1 and 1; diagonal values vary lie in the range of 0.6 to 1.4 times the base value
    // off diagonal values vary between 0.7 to 1.1 of the corresponding diagonal value (symmetrization is conducted)
    m_cont_freq_matrix_variable =
        std::move(std::make_unique<ContactFrequencyVariableElement>(ContactFrequencyVariableElement{
            cont_freq_matrix,
            nb_damp_dist,
            day_dist,
            damp_diag_base_dist,
            damp_diag_rel_dist,
            damp_offdiag_rel_dist}));
	return m_cont_freq_matrix_variable;
}
	

void write_parameters(std::vector<epi::SecirParams> const& params, epi::ContactFrequencyMatrix const& cont_freq_matrix, double t0, double tmax, double dt, int runs, const dist_params& dists, const std::string& filename)
{
	char* docString = NULL;
	double* myFloat = NULL;
	
	TixiDocumentHandle handle;




	const int nb_groups = params.size();
	const int nb_damp = cont_freq_matrix.get_dampings(0, 0).get_dampings_vector().size();
	//create xml doc
	tixiCreateDocument("Parameters", &handle);
	tixiAddIntegerElement(handle, "/Parameters", "Runs", runs, "%d");
	tixiAddIntegerElement(handle, "/Parameters", "NumberOfGroups", nb_groups, "%d");
	tixiAddDoubleElement(handle, "/Parameters", "T0", t0, "%g");
	tixiAddDoubleElement(handle, "/Parameters", "TMax", tmax, "%g");
	tixiAddDoubleElement(handle, "/Parameters", "dt", dt, "%g");
	tixiAddIntegerElement(handle, "/Parameters", "NumberOfDampings", nb_damp, "%d");

	//add float vector
	double** contact = new double*[nb_groups];
	for (int i = 0; i < nb_groups; ++i)
	{
		contact[i] = new double[nb_groups];
		for (int j = 0; j < nb_groups; ++j)
		{
			contact[i][j] = cont_freq_matrix.get_cont_freq(i, j);
		}
		std::string name = "ContactRateGroup" + std::to_string(i + 1);
		tixiAddFloatVector(handle, "/Parameters", name.c_str(), contact[i], nb_groups, "%g");
	}


	tixiCreateElement(handle, "/Parameters", "Dampings");
	double damping;
	for (int i = 0; i < nb_damp;++i)
	{
		tixiCreateElement(handle, "/Parameters/Dampings", ("Damp" + std::to_string(i + 1)).c_str());
		tixiAddDoubleElement(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1)).c_str(), "Day", cont_freq_matrix.get_dampings(0, 0).get_dampings_vector().at(i).day, "%g");
		for (int k = 0; k < nb_groups;++k)
		{
			damping = cont_freq_matrix.get_dampings(k, k).get_dampings_vector().at(i).factor;
			std::string name = "Group" + std::to_string(k + 1);
			tixiAddDoubleElement(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1)).c_str(), name.c_str(), damping, "%g");
		}
	}



	for (int i = 0; i < nb_groups; ++i)
	{

		std::string name = "Group" + std::to_string(i + 1);
		std::string path = "/Parameters/" + name;
		std::string pop_path = path + "/Population";
		std::string time_path = path + "/StageTimes";
		std::string prob_path = path + "/Probabilities";

		tixiCreateElement(handle, "/Parameters", name.c_str());
		tixiCreateElement(handle, path.c_str(), "Population");
		tixiCreateElement(handle, path.c_str(), "StageTimes");
		tixiCreateElement(handle, path.c_str(), "Probabilities");


		write_dist(handle, pop_path, "Total", dists.dist_total, params[i].populations.get_total_t0(), dists.total);
		write_dist(handle, pop_path, "Exposed", dists.dist_exposed, params[i].populations.get_exposed_t0(), dists.exposed);
		write_dist(handle, pop_path, "Carrier", dists.dist_carrier, params[i].populations.get_carrier_t0(), dists.carrier);
		write_dist(handle, pop_path, "Infectious", dists.dist_infectious, params[i].populations.get_infectious_t0(), dists.infectious);
		write_dist(handle, pop_path, "Hospitalized", dists.dist_hospital, params[i].populations.get_hospitalized_t0(), dists.hospital);
		write_dist(handle, pop_path, "ICU", dists.dist_icu, params[i].populations.get_icu_t0(), dists.icu);
		write_dist(handle, pop_path, "Recovered", dists.dist_recovered, params[i].populations.get_recovered_t0(), dists.recovered);
		write_dist(handle, pop_path, "Dead", dists.dist_dead, params[i].populations.get_dead_t0(), dists.dead);
		
		write_dist(handle, time_path, "Incubation", dists.dist_tinc, 1.0 / params[i].times.get_incubation_inv(), dists.tinc);
		write_dist(handle, time_path, "InfectiousMild", dists.dist_tinfmild, 1.0 / params[i].times.get_infectious_mild_inv(), dists.tinfmild);
		write_dist(handle, time_path, "SerialInterval", dists.dist_tserint, 1.0 / params[i].times.get_serialinterval_inv();, dists.tserint);
		write_dist(handle, time_path, "HospitalizedToHome", dists.dist_thosp2home, 1.0 / params[i].times.get_hospitalized_to_home_inv(), dists.thosp2home);
		write_dist(handle, time_path, "HomeToHospitalized", dists.dist_thome2hosp, 1.0 / params[i].times.get_home_to_hospitalized_inv(), dists.thome2hosp);
		write_dist(handle, time_path, "HospitalizedToICU", dists.dist_thosp2icu, 1.0 / params[i].times.get_hospitalized_to_icu_inv(), dists.thosp2icu);
		write_dist(handle, time_path, "ICUToHome", dists.dist_ticu2home, 1.0 / params[i].times.get_icu_to_home_inv(), dists.ticu2home);
		write_dist(handle, time_path, "InfectiousAsymp", dists.dist_tinfasy, 1.0 / params[i].times.get_infectious_asymp_inv(), dists.tinfasy);
		write_dist(handle, time_path, "ICUToDeath", dists.dist_ticu2death, 1.0 / params[i].times.get_icu_to_dead_inv(), dists.ticu2death);
		
		write_dist(handle, prob_path, "InfectionFromContact", dists.inf_cont, params[i].probabilities.get_infection_from_contact(), dists.inf_cont);
		write_dist(handle, prob_path, "AsympPerInfection", dists.dist_alpha, params[i].probabilities.get_asymp_per_infectious(), dists.alpha);
		write_dist(handle, prob_path, "RiskFromSymptomatic", dists.dist_beta, params[i].probabilities.get_risk_from_symptomatic(), dists.beta);
		write_dist(handle, prob_path, "HospitalizedPerInfectious", dists.dist_rho, params[i].probabilities.get_hospitalized_per_infectious(), dists.rho);
		write_dist(handle, prob_path, "ICUPerHospitalized", dists.dist_theta, params[i].probabilities.get_icu_per_hospitalized(), dists.theta);
		write_dist(handle, prob_path, "DeadPerICU", dists.dist_delta, params[i].probabilities.get_dead_per_icu(), dists.delta);

		
	}

	delete contact;


	//write to file
	tixiSaveDocument(handle, filename.c_str());
	tixiCloseDocument(handle);
	handle = 0;

	tixiCleanup();

}

double draw_number(double* values, char* dist)
{
    if (strcmp("uniform", dist) == 0)
	{
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(values[1], values[2]);
		double number = distribution(generator);
		return number;
	}
	else if (strcmp("normal", dist) == 0)
	{
		std::random_device mch;
		std::default_random_engine generator(mch());
		std::normal_distribution<double> distribution(values[0], values[3]);
		double number = distribution(generator);
		if (number > values[2])
		{
			number = values[2];
		}
		else if (number < values[1])
		{
			number = values[1];
		}
		return number;
	}
    else
    {
        return values[0];
    }
}




file read_parameters(const std::string& filename)
{
	TixiDocumentHandle handle = -1;
	tixiOpenDocument(filename.c_str(), &handle);
	int nb_groups;
	int nb_damp;
	int runs;
	double t0;
	double tmax;
	double dt;

	char* dist;
	int dist_size;

	tixiGetIntegerElement(handle, "/Parameters/NumberOfGroups", &nb_groups);
	tixiGetIntegerElement(handle, "/Parameters/NumberOfDampings", &nb_damp);
	tixiGetIntegerElement(handle, "/Parameters/Runs", &runs);
	tixiGetTextElement(handle, "/Parameters/Distribution", &dist);
	tixiGetDoubleElement(handle, "/Parameters/T0", &t0);
	tixiGetDoubleElement(handle, "/Parameters/TMax", &tmax);
	tixiGetDoubleElement(handle, "/Parameters/dt", &dt);

	

	if (strcmp("none", dist) == 0)
	{
		dist_size = 0;
	}
	else if (strcmp("normal", dist) == 0)
	{
		dist_size = 3;
	}
	else if (strcmp("uniform", dist) == 0)
	{
		dist_size = 2;
	}

	std::vector<std::vector<epi::SecirParams>> all_params;
	std::vector<epi::ContactFrequencyMatrix> all_contact_freq_matrix;

	for (int i = 0; i < runs; ++i)
	{



		std::vector<epi::SecirParams> params{ epi::SecirParams{} };
		epi::ContactFrequencyMatrix contact_freq_matrix{ (size_t)nb_groups };
		for (size_t i = 1; i < nb_groups; i++) {
			params.push_back(epi::SecirParams{});
		}

		for (size_t i = 0; i < nb_groups; i++) {

			std::string name = "Group" + std::to_string(i + 1);
			std::string path = "/Parameters/" + name;
			std::string pop_path = path + "/Population";
			std::string time_path = path + "/StageTimes";
			std::string prob_path = path + "/Probabilities";

			double* contact = new double[nb_groups];
			tixiGetFloatVector(handle, ("/Parameters/ContactRateGroup" + std::to_string(i + 1)).c_str(), &contact, nb_groups);

			for (int j = 0; j < nb_groups; ++j)
			{
				contact_freq_matrix.set_cont_freq(contact[j], i, j);
			}


			double* tinc = new double[dist_size + 1];
			double* tinfmild = new double[dist_size + 1];
			double* tserint = new double[dist_size + 1];
			double* thosp2home = new double[dist_size + 1];
			double* thome2hosp = new double[dist_size + 1];
			double* thosp2icu = new double[dist_size + 1];
			double* ticu2home = new double[dist_size + 1];
			double* tinfasy = new double[dist_size + 1];
			double* ticu2death = new double[dist_size + 1];

			double* inf_cont = new double[dist_size + 1];
			double* alpha = new double[dist_size + 1];
			double* beta = new double[dist_size + 1];
			double* rho = new double[dist_size + 1];
			double* theta = new double[dist_size + 1];
			double* delta = new double[dist_size + 1];

			double nb_total_t0;
			double nb_exp_t0;
			double nb_car_t0;
			double nb_inf_t0;
			double nb_hosp_t0;
			double nb_icu_t0;
			double nb_rec_t0;
			double nb_dead_t0;

			



			tixiGetFloatVector(handle, (time_path + "/Incubation").c_str(), &tinc, dist_size + 1);
			tixiGetFloatVector(handle, (time_path + "/InfectiousMild").c_str(), &tinfmild, dist_size + 1);
			tixiGetFloatVector(handle, (time_path + "/SerialInterval").c_str(), &tserint, dist_size + 1);
			tixiGetFloatVector(handle, (time_path + "/HospitalizedToHome").c_str(), &thosp2home, dist_size + 1);
			tixiGetFloatVector(handle, (time_path + "/HomeToHospitalized").c_str(), &thome2hosp, dist_size + 1);
			tixiGetFloatVector(handle, (time_path + "/HospitalizedToICU").c_str(), &thosp2icu, dist_size + 1);
			tixiGetFloatVector(handle, (time_path + "/ICUToHome").c_str(), &ticu2home, dist_size + 1);
			tixiGetFloatVector(handle, (time_path + "/InfectiousAsymp").c_str(), &tinfasy, dist_size + 1);
			tixiGetFloatVector(handle, (time_path + "/ICUToDeath").c_str(), &ticu2death, dist_size + 1);

			tixiGetDoubleElement(handle, (pop_path + "/Total").c_str(), &nb_total_t0);
			tixiGetDoubleElement(handle, (pop_path + "/Exposed").c_str(), &nb_exp_t0);
			tixiGetDoubleElement(handle, (pop_path + "/Carrier").c_str(), &nb_car_t0);
			tixiGetDoubleElement(handle, (pop_path + "/Infectious").c_str(), &nb_inf_t0);
			tixiGetDoubleElement(handle, (pop_path + "/Hospital").c_str(), &nb_hosp_t0);
			tixiGetDoubleElement(handle, (pop_path + "/ICU").c_str(), &nb_icu_t0);
			tixiGetDoubleElement(handle, (pop_path + "/Recovered").c_str(), &nb_rec_t0);
			tixiGetDoubleElement(handle, (pop_path + "/Dead").c_str(), &nb_dead_t0);

			tixiGetFloatVector(handle, (prob_path + "/InfectionFromContact").c_str(), &inf_cont, dist_size + 1);
			tixiGetFloatVector(handle, (prob_path + "/AsympPerInfection").c_str(), &alpha, dist_size + 1);
			tixiGetFloatVector(handle, (prob_path + "/RiskFromSymptomatic").c_str(), &beta, dist_size + 1);
			tixiGetFloatVector(handle, (prob_path + "/HospitalizedPerInfectious").c_str(), &rho, dist_size + 1);
			tixiGetFloatVector(handle, (prob_path + "/ICUPerHospitalized").c_str(), &theta, dist_size + 1);
			tixiGetFloatVector(handle, (prob_path + "/DeadPerICU").c_str(), &delta, dist_size + 1);



			params[i].times.set_incubation(draw_number(tinc, dist));
			params[i].times.set_infectious_mild(draw_number(tinfmild, dist));
			params[i].times.set_serialinterval(draw_number(tserint, dist));
			params[i].times.set_hospitalized_to_home(draw_number(thosp2home, dist));
			params[i].times.set_home_to_hospitalized(draw_number(thome2hosp, dist));
			params[i].times.set_hospitalized_to_icu(draw_number(thosp2icu, dist));
			params[i].times.set_icu_to_home(draw_number(ticu2home, dist));
			params[i].times.set_infectious_asymp(draw_number(tinfasy, dist));
			params[i].times.set_icu_to_death(draw_number(ticu2death, dist));

			params[i].populations.set_total_t0(nb_total_t0);
			params[i].populations.set_exposed_t0(nb_exp_t0);
			params[i].populations.set_carrier_t0(nb_car_t0);
			params[i].populations.set_infectious_t0(nb_inf_t0);
			params[i].populations.set_hospital_t0(nb_hosp_t0);
			params[i].populations.set_icu_t0(nb_icu_t0);
			params[i].populations.set_recovered_t0(nb_rec_t0);
			params[i].populations.set_dead_t0(nb_dead_t0);

			params[i].probabilities.set_infection_from_contact(draw_number(inf_cont, dist));
			params[i].probabilities.set_asymp_per_infectious(draw_number(alpha, dist));
			params[i].probabilities.set_risk_from_symptomatic(draw_number(beta, dist));
			params[i].probabilities.set_hospitalized_per_infectious(draw_number(rho, dist));
			params[i].probabilities.set_icu_per_hospitalized(draw_number(theta, dist));
			params[i].probabilities.set_dead_per_icu(draw_number(delta, dist));
		}

		for (int i = 0; i < nb_damp; ++i)
		{
			double day;
			double* damping = new double[nb_groups];

			tixiGetDoubleElement(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1) + "/Day").c_str(), &day);
			for (int k = 0; k < nb_groups; ++k)
			{
				tixiGetFloatVector(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1) + "/Group" + std::to_string(k + 1)).c_str(), &damping, nb_groups);
				for (int l = 0; l < nb_groups; ++l)
				{
					epi::Damping dummy(day, damping[l]);
					contact_freq_matrix.add_damping(dummy, k, l);
				}
			}

		}

		all_params.push_back(params);
		all_contact_freq_matrix.push_back(contact_freq_matrix);
	}
	tixiCloseDocument(handle);

	return { nb_groups, runs, t0, tmax, dt, all_params, all_contact_freq_matrix };
}


parameter_space_t::parameter_space_t(std::string& parameter_filename)
{
	TixiDocumentHandle handle = -1;
	tixiOpenDocument(filename.c_str(), &handle);
	int nb_groups;
	int nb_damp;
	int runs;
	double t0;
	double tmax;
	double dt;


	tixiGetIntegerElement(handle, "/Parameters/NumberOfGroups", &nb_groups);
	tixiGetIntegerElement(handle, "/Parameters/NumberOfDampings", &nb_damp);
	tixiGetIntegerElement(handle, "/Parameters/Runs", &runs);
	tixiGetDoubleElement(handle, "/Parameters/T0", &t0);
	tixiGetDoubleElement(handle, "/Parameters/TMax", &tmax);
	tixiGetDoubleElement(handle, "/Parameters/dt", &dt);
	
	epi::ContactFrequencyMatrix contact_freq_matrix{ (size_t)nb_groups };

	
	
	for (size_t i = 0; i < nb_groups; i++) {

		std::string name = "Group" + std::to_string(i + 1);
		std::string path = "/Parameters/" + name;
		std::string pop_path = path + "/Population";
		std::string time_path = path + "/StageTimes";
		std::string prob_path = path + "/Probabilities";

		double* contact = new double[nb_groups];
		tixiGetFloatVector(handle, ("/Parameters/ContactRateGroup" + std::to_string(i + 1)).c_str(), &contact, nb_groups);

		for (int j = 0; j < nb_groups; ++j)
		{
			contact_freq_matrix.set_cont_freq(contact[j], i, j);
		}

        // fixed size groups
        // total
		
		double total_t0;
		double hosp_t0;
		double icu_t0;
		double dead_t0;
		
		tixiGetDoubleElement(handle, (pop_path + "/Total/Mean").c_str(), &total_t0);
		tixiGetDoubleElement(handle, (pop_path + "/Hospitalized/Mean").c_str(), &hosp_t0);
		tixiGetDoubleElement(handle, (pop_path + "/ICU/Mean").c_str(), &icu_t0);
		tixiGetDoubleElement(handle, (pop_path + "/Dead/Mean").c_str(), &dead_t0);
		
        m_total.push_back(total_t0);
        m_hospitalized.push_back(hosp_t0);
        m_icu.push_back(icu_t0);
        m_dead.push_back(dead_t0);

        // variably sized groups
        // exposed
        m_exposed.push_back(read_dist(handle, pop_path + "/Exposed"));

        // carrier
        value_params = params[i].populations.get_carrier_t0();
        m_carrier.push_back(read_dist(handle, pop_path + "/Carrier"));

        // infectious
        m_infectious.push_back(read_dist(handle, pop_path + "/Infectious"));

        // recovered
        m_recovered.push_back(read_dist(handle, pop_path + "/Recovered"));
		
		
		// times
		// incubation time
        m_incubation.push_back(read_dist(handle, time_path + "/Incubation"));

        // infectious time (mild)
        m_inf_mild.push_back(read_dist(handle, time_path + "/InfectiousMild"));

        // serial interval
        m_serial_int.push_back(read_dist(handle, time_path + "/SerialInterval"));

        // infectious to recovered (hospitalized to home)
        m_hosp_to_rec.push_back(read_dist(handle, time_path + "/HospitalizedToHome"));

        // infectious to hospitalized (home to hospitalized)
        m_inf_to_hosp.push_back(read_dist(handle, time_path + "/HomeToHospitalized"));

        // infectious (asymptomatic)
        m_inf_asymp.push_back(read_dist(handle, time_path + "/InfectiousAsymp"));

        // hospitalized to ICU
        m_hosp_to_icu.push_back(read_dist(handle, time_path + "/HospitalizedToICU"));

        // ICU to recovered
        m_icu_to_rec.push_back(read_dist(handle, time_path + "/ICUToHome"));

        // ICU to death
        m_icu_to_death.push_back(read_dist(handle, time_path + "/ICUToDeath"));
		
		
		// probabilities
		// infection from contact
        m_inf_from_cont.push_back(read_dist(handle, prob_path + "/InfectionFromContact"));

        // asymptomatic per infectious
        m_asymp_per_inf.push_back(read_dist(handle, prob_path + "/AsympPerInfection"));

        // risk of infection from infectious
        m_risk_from_symp.push_back(read_dist(handle, prob_path + "/RiskFromSymptomatic"));

        // deaths per icu treatments
        m_death_per_icu.push_back(read_dist(handle, prob_path + "/DeadPerICU"));

        // hospitalized per infections
        m_hosp_per_inf.push_back(read_dist(handle, prob_path + "/HospitalizedPerInfectious"));

        // icu treatments per hospitalized
        m_icu_per_hosp.push_back(read_dist(handle, prob_path + "/ICUPerHospitalized"));
		
		for (int j = 0; j < nb_groups; ++j){
			contact_freq_matrix.set_cont_freq(contact[j], i, j);
		}
		
		
    }
	

	 nb_damp_dist = std::make_unique<ParameterDistributionUniform>(ParameterDistributionUniform(1, (tmax - t0) / 10));
	 day_dist = std::make_unique<ParameterDistributionUniform>(ParameterDistributionUniform(t0, tmax));
	 damp_diag_base_dist = std::make_unique<ParameterDistributionUniform>(ParameterDistributionUniform(0.1, 1));
	 damp_diag_rel_dist = std::make_unique<ParameterDistributionUniform>(ParameterDistributionUniform(0.6, 1.4));
	 damp_offdiag_rel_dist = std::make_unique<ParameterDistributionUniform>(ParameterDistributionUniform(0.7, 1.1));
	 
	 nb_damp_dist.add_predefined_sample(nb_damp);
	 for (int i = 0; i < nb_damp; ++i)
	 {
		double day;
		double damping;

		tixiGetDoubleElement(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1) + "/Day").c_str(), &day);
		day_dist.add_predefined_sample(day);
		for (int k = 0; k < nb_groups; ++k)
		{
			tixiGetDoubleElement(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1) + "/Group" + std::to_string(k + 1)).c_str(), &damping);
			damp_diag_base.add_predefined_sample(damping);
			damp_diag_rel_dist.add_predefined_sample(1.0);
			for (int l = k+1; l<nb_groups; ++l)
				damp_offdiag_rel_dist.add_predefined_sample(1.0);
		}
	 }
	 
	 
	// maximum number of dampings; to avoid overfitting only allow one damping for every 10 days simulated
    // damping base values are between 0.1 and 1; diagonal values vary lie in the range of 0.6 to 1.4 times the base value
    // off diagonal values vary between 0.7 to 1.1 of the corresponding diagonal value (symmetrization is conducted)
    m_cont_freq_matrix_variable =
        std::move(std::make_unique<ContactFrequencyVariableElement>(ContactFrequencyVariableElement{
            cont_freq_matrix,
            nb_damp_dist,
            day_dist,
            damp_diag_base_dist,
            damp_diag_rel_dist,
            damp_offdiag_rel_dist}));
    // TODO: implement
    assert(0 && "This function not implemented yet and needs a file read method.");
}


}
