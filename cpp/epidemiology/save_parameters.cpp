#include <epidemiology/secir.h>
#include <epidemiology/damping.h>
#include <vector>

#include <iostream>
#include <string>

#include <tixi.h>

void write_parameters(std::vector<epi::SecirParams> const& params, epi::ContactFrequencyMatrix const& cont_freq_matrix, double t0, double tmax, double dt, std::string filename)
{
	char* docString = NULL;
	double* myFloat = NULL;
	int i = 0;
	TixiDocumentHandle handle;

	printf("TIXI Demo\n\n");
	printf("Using TIXI %s\n", tixiGetVersion());

	const int nb_groups = params.size();
	const int nb_damp = cont_freq_matrix.get_dampings(0, 0).get_dampings_vector().size();
	//create xml doc
	tixiCreateDocument("Parameters", &handle);
	//tixiAddHeader(handle, "tixiDemo", tixiGetVersion(), "Martin Siggel");
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

		tixiAddDoubleElement(handle, time_path.c_str(), "Incubation", 1.0 / params[i].times.get_incubation_inv(), "%g");
		tixiAddDoubleElement(handle, time_path.c_str(), "InfectiousMild", 1.0 / params[i].times.get_infectious_mild_inv(), "%g");
		tixiAddDoubleElement(handle, time_path.c_str(), "SerialInterval", 1.0 / params[i].times.get_serialinterval_inv(), "%g");
		tixiAddDoubleElement(handle, time_path.c_str(), "HospitalizedToHome", 1.0 / params[i].times.get_hospitalized_to_home_inv(), "%g");
		tixiAddDoubleElement(handle, time_path.c_str(), "HomeToHospitalized", 1.0 / params[i].times.get_home_to_hospitalized_inv(), "%g");
		tixiAddDoubleElement(handle, time_path.c_str(), "HospitalizedToICU", 1.0 / params[i].times.get_hospitalized_to_icu_inv(), "%g");
		tixiAddDoubleElement(handle, time_path.c_str(), "ICUToHome", 1.0 / params[i].times.get_icu_to_home_inv(), "%g");
		tixiAddDoubleElement(handle, time_path.c_str(), "InfectiousAsymp", 1.0 / params[i].times.get_infectious_asymp_inv(), "%g");
		tixiAddDoubleElement(handle, time_path.c_str(), "ICUToDeath", 1.0 / params[i].times.get_icu_to_dead_inv(), "%g");

		tixiAddDoubleElement(handle, prob_path.c_str(), "InfectionFromContact", params[i].probabilities.get_infection_from_contact(), "%g");
		tixiAddDoubleElement(handle, prob_path.c_str(), "AsympPerInfection", params[i].probabilities.get_asymp_per_infectious(), "%g");
		tixiAddDoubleElement(handle, prob_path.c_str(), "RiskFromSymptomatic", params[i].probabilities.get_risk_from_symptomatic(), "%g");
		tixiAddDoubleElement(handle, prob_path.c_str(), "HospitalizedPerInfectious", params[i].probabilities.get_hospitalized_per_infectious(), "%g");
		tixiAddDoubleElement(handle, prob_path.c_str(), "ICUPerHospitalized", params[i].probabilities.get_icu_per_hospitalized(), "%g");
		tixiAddDoubleElement(handle, prob_path.c_str(), "DeadPerICU", params[i].probabilities.get_dead_per_icu(), "%g");
	}

	delete contact;


	//write to file
	tixiSaveDocument(handle, filename.c_str());
	tixiCloseDocument(handle);
	handle = 0;

	tixiCleanup();

}

struct file
{
	int nb_groups;
	double t0;
	double tmax;
	double dt;
	std::vector<epi::SecirParams> params;
	epi::ContactFrequencyMatrix contact_freq_matrix;
};

file read_parameters(std::string filename)
{
	TixiDocumentHandle handle = -1;
	tixiOpenDocument(filename.c_str(), &handle);
	int nb_groups;
	int nb_damp;
	double t0;
	double tmax;
	double dt;

	tixiGetIntegerElement(handle, "/Parameters/NumberOfGroups", &nb_groups);
	tixiGetIntegerElement(handle, "/Parameters/NumberOfDampings", &nb_damp);
	tixiGetDoubleElement(handle, "/Parameters/T0", &t0);
	tixiGetDoubleElement(handle, "/Parameters/TMax", &tmax);
	tixiGetDoubleElement(handle, "/Parameters/dt", &dt);

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


		double tinc;
		double tinfmild;
		double tserint;
		double thosp2home;
		double thome2hosp;
		double thosp2icu;
		double ticu2home;
		double tinfasy;
		double ticu2death;

		double nb_total_t0;
		double nb_exp_t0;
		double nb_car_t0;
		double nb_inf_t0;
		double nb_hosp_t0;
		double nb_icu_t0;
		double nb_rec_t0;
		double nb_dead_t0;

		double inf_cont;
		double alpha;
		double beta;
		double rho;
		double theta;
		double delta;



		tixiGetDoubleElement(handle, (time_path + "/Incubation").c_str(), &tinc);
		tixiGetDoubleElement(handle, (time_path + "/InfectiousMild").c_str(), &tinfmild);
		tixiGetDoubleElement(handle, (time_path + "/SerialInterval").c_str(), &tserint);
		tixiGetDoubleElement(handle, (time_path + "/HospitalizedToHome").c_str(), &thosp2home);
		tixiGetDoubleElement(handle, (time_path + "/HomeToHospitalized").c_str(), &thome2hosp);
		tixiGetDoubleElement(handle, (time_path + "/HospitalizedToICU").c_str(), &thosp2icu);
		tixiGetDoubleElement(handle, (time_path + "/ICUToHome").c_str(), &ticu2home);
		tixiGetDoubleElement(handle, (time_path + "/InfectiousAsymp").c_str(), &tinfasy);
		tixiGetDoubleElement(handle, (time_path + "/ICUToDeath").c_str(), &ticu2death);

		tixiGetDoubleElement(handle, (pop_path + "/Total").c_str(), &nb_total_t0);
		tixiGetDoubleElement(handle, (pop_path + "/Exposed").c_str(), &nb_exp_t0);
		tixiGetDoubleElement(handle, (pop_path + "/Carrier").c_str(), &nb_car_t0);
		tixiGetDoubleElement(handle, (pop_path + "/Infectious").c_str(), &nb_inf_t0);
		tixiGetDoubleElement(handle, (pop_path + "/Hospital").c_str(), &nb_hosp_t0);
		tixiGetDoubleElement(handle, (pop_path + "/ICU").c_str(), &nb_icu_t0);
		tixiGetDoubleElement(handle, (pop_path + "/Recovered").c_str(), &nb_rec_t0);
		tixiGetDoubleElement(handle, (pop_path + "/Dead").c_str(), &nb_dead_t0);

		tixiGetDoubleElement(handle, (prob_path + "/InfectionFromContact").c_str(), &inf_cont);
		tixiGetDoubleElement(handle, (prob_path + "/AsympPerInfection").c_str(), &alpha);
		tixiGetDoubleElement(handle, (prob_path + "/RiskFromSymptomatic").c_str(), &beta);
		tixiGetDoubleElement(handle, (prob_path + "/HospitalizedPerInfectious").c_str(), &rho);
		tixiGetDoubleElement(handle, (prob_path + "/ICUPerHospitalized").c_str(), &theta);
		tixiGetDoubleElement(handle, (prob_path + "/DeadPerICU").c_str(), &delta);




		params[i].times.set_incubation(tinc);
		params[i].times.set_infectious_mild(tinfmild);
		params[i].times.set_serialinterval(tserint);
		params[i].times.set_hospitalized_to_home(thosp2home);
		params[i].times.set_home_to_hospitalized(thome2hosp);
		params[i].times.set_hospitalized_to_icu(thosp2icu);
		params[i].times.set_icu_to_home(ticu2home);
		params[i].times.set_infectious_asymp(tinfasy);
		params[i].times.set_icu_to_death(ticu2death);

		params[i].populations.set_total_t0(nb_total_t0);
		params[i].populations.set_exposed_t0(nb_exp_t0);
		params[i].populations.set_carrier_t0(nb_car_t0);
		params[i].populations.set_infectious_t0(nb_inf_t0);
		params[i].populations.set_hospital_t0(nb_hosp_t0);
		params[i].populations.set_icu_t0(nb_icu_t0);
		params[i].populations.set_recovered_t0(nb_rec_t0);
		params[i].populations.set_dead_t0(nb_dead_t0);

		params[i].probabilities.set_infection_from_contact(1.0);
		params[i].probabilities.set_asymp_per_infectious(alpha);
		params[i].probabilities.set_risk_from_symptomatic(beta);
		params[i].probabilities.set_hospitalized_per_infectious(rho);
		params[i].probabilities.set_icu_per_hospitalized(theta);
		params[i].probabilities.set_dead_per_icu(delta);
	}

	for (int i = 0; i < nb_damp; ++i)
	{
		double day;
		double* damping = new double[nb_groups];

		tixiGetDoubleElement(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1) + "/Day").c_str(), &day);
		for (int k = 0; k < nb_groups; ++k)
		{
			tixiGetFloatVector(handle, ("/Parameters/Dampings/Damp" + std::to_string(i + 1) + "/Group" + std::to_string(k+1)).c_str(), &damping, nb_groups);
			for (int l = 0; l < nb_groups; ++l)
			{
				epi::Damping dummy(day, damping[l]);
				contact_freq_matrix.add_damping(dummy, k, l);
			}
		}
		
	}



	tixiCloseDocument(handle);

	return { nb_groups, t0, tmax, dt, params, contact_freq_matrix };
}
