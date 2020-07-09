#include <epidemiology/secir.h>
#include <epidemiology/damping.h>
#include <vector>

#include <iostream>
#include <string>
#include <random>

#include <tixi.h>

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

void write_parameters(std::vector<epi::SecirParams> const& params, epi::ContactFrequencyMatrix const& cont_freq_matrix, double t0, double tmax, double dt, int runs, const std::string& dist, const dist_params& dists, const std::string& filename)
{
	char* docString = NULL;
	double* myFloat = NULL;
	int i = 0;
	int dist_size;
	if (dist == "none")
	{
		dist_size = 0;
	}
	else if (dist == "normal")
	{
		dist_size = 3;
	}
	else if (dist == "uniform")
	{
		dist_size = 2;
	}
	else
	{
		dist_size = 0;
	}
	
	TixiDocumentHandle handle;




	const int nb_groups = params.size();
	const int nb_damp = cont_freq_matrix.get_dampings(0, 0).get_dampings_vector().size();
	//create xml doc
	tixiCreateDocument("Parameters", &handle);
	tixiAddIntegerElement(handle, "/Parameters", "Runs", runs, "%d");
	tixiAddTextElement(handle, "/Parameters", "Distribution", dist.c_str());
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
	for (int i = 0; i < nb_damp;++i)
	{
		tixiCreateElement(handle, "/Parameters/Dampings", ("Damp" + std::to_string(i + 1)).c_str());
		tixiAddDoubleElement(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1)).c_str(), "Day", cont_freq_matrix.get_dampings(0, 0).get_dampings_vector().at(i).day, "%g");
		for (int k = 0; k < nb_groups;++k)
		{
			for (int l = 0; l < nb_groups; ++l)
			{
				contact[k][l] = cont_freq_matrix.get_dampings(k, l).get_dampings_vector().at(i).factor;
			}
			std::string name = "Group" + std::to_string(k + 1);
			tixiAddFloatVector(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1)).c_str(), name.c_str(), contact[k], nb_groups, "%g");
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



		tixiAddDoubleElement(handle, pop_path.c_str(), "Total", params[i].populations.get_total_t0(), "%g");
		tixiAddDoubleElement(handle, pop_path.c_str(), "Exposed", params[i].populations.get_exposed_t0(), "%g");
		tixiAddDoubleElement(handle, pop_path.c_str(), "Carrier", params[i].populations.get_carrier_t0(), "%g");
		tixiAddDoubleElement(handle, pop_path.c_str(), "Infectious", params[i].populations.get_infectious_t0(), "%g");
		tixiAddDoubleElement(handle, pop_path.c_str(), "Hospital", params[i].populations.get_hospitalized_t0(), "%g");
		tixiAddDoubleElement(handle, pop_path.c_str(), "ICU", params[i].populations.get_icu_t0(), "%g");
		tixiAddDoubleElement(handle, pop_path.c_str(), "Recovered", params[i].populations.get_recovered_t0(), "%g");
		tixiAddDoubleElement(handle, pop_path.c_str(), "Dead", params[i].populations.get_dead_t0(), "%g");


		double* tinc = new double[dist_size +1];
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

		tinc[0] = 1.0 / params[i].times.get_incubation_inv();
		tinfmild[0] = 1.0 / params[i].times.get_infectious_mild_inv();
		tserint[0] = 1.0 / params[i].times.get_serialinterval_inv();
		thosp2home[0] = 1.0 / params[i].times.get_hospitalized_to_home_inv();
		thome2hosp[0] = 1.0 / params[i].times.get_home_to_hospitalized_inv();
		thosp2icu[0] = 1.0 / params[i].times.get_hospitalized_to_icu_inv();
		ticu2home[0] = 1.0 / params[i].times.get_icu_to_home_inv();
		tinfasy[0] = 1.0 / params[i].times.get_infectious_asymp_inv();
		ticu2death[0] = 1.0 / params[i].times.get_icu_to_dead_inv();

		inf_cont[0] = params[i].probabilities.get_infection_from_contact();
		alpha[0] = params[i].probabilities.get_asymp_per_infectious();
		beta[0] = params[i].probabilities.get_risk_from_symptomatic();
		rho[0] = params[i].probabilities.get_hospitalized_per_infectious();
		theta[0] = params[i].probabilities.get_icu_per_hospitalized();
		delta[0] = params[i].probabilities.get_dead_per_icu();
		for (int i = 0; i < dist_size; ++i)
		{
			tinc[i + 1] = dists.tinc[i];
			tinfmild[i + 1] = dists.tinfmild[i];
			tserint[i + 1] = dists.tserint[i];
			thosp2home[i + 1] = dists.thosp2home[i];
			thome2hosp[i + 1] = dists.thome2hosp[i];
			thosp2icu[i + 1] = dists.thosp2icu[i];
			ticu2home[i + 1] = dists.ticu2home[i];
			tinfasy[i + 1] = dists.tinfasy[i];
			ticu2death[i + 1] = dists.ticu2death[i];

			inf_cont[i + 1] = dists.inf_cont[i];
			alpha[i + 1] = dists.alpha[i];
			beta[i + 1] = dists.beta[i];
			rho[i + 1] = dists.rho[i];
			theta[i + 1] = dists.theta[i];
			delta[i + 1] = dists.delta[i];
		}
		


		tixiAddFloatVector(handle, time_path.c_str(), "Incubation", tinc, dist_size+1, "%g");
		tixiAddFloatVector(handle, time_path.c_str(), "InfectiousMild", tinfmild, dist_size + 1, "%g");
		tixiAddFloatVector(handle, time_path.c_str(), "SerialInterval",tserint, dist_size + 1, "%g");
		tixiAddFloatVector(handle, time_path.c_str(), "HospitalizedToHome", thosp2home, dist_size + 1, "%g");
		tixiAddFloatVector(handle, time_path.c_str(), "HomeToHospitalized", thome2hosp, dist_size + 1, "%g");
		tixiAddFloatVector(handle, time_path.c_str(), "HospitalizedToICU", thosp2icu, dist_size + 1, "%g");
		tixiAddFloatVector(handle, time_path.c_str(), "ICUToHome", ticu2home, dist_size + 1, "%g");
		tixiAddFloatVector(handle, time_path.c_str(), "InfectiousAsymp", tinfasy, dist_size + 1, "%g");
		tixiAddFloatVector(handle, time_path.c_str(), "ICUToDeath", ticu2death, dist_size + 1, "%g");

		tixiAddFloatVector(handle, prob_path.c_str(), "InfectionFromContact", inf_cont, dist_size + 1, "%g");
		tixiAddFloatVector(handle, prob_path.c_str(), "AsympPerInfection", alpha, dist_size + 1, "%g");
		tixiAddFloatVector(handle, prob_path.c_str(), "RiskFromSymptomatic", beta, dist_size + 1, "%g");
		tixiAddFloatVector(handle, prob_path.c_str(), "HospitalizedPerInfectious", rho, dist_size + 1, "%g");
		tixiAddFloatVector(handle, prob_path.c_str(), "ICUPerHospitalized", theta, dist_size + 1, "%g");
		tixiAddFloatVector(handle, prob_path.c_str(), "DeadPerICU", delta, dist_size + 1, "%g");
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

}
