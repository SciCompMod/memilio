#ifndef SEIR_PARAM
#define SEIR_PARAM

/**
 * This defined a damping factor for a
 * mitigation strategy for one point in time.
 */
template <typename T>
struct damping {
    T day;
    T factor;

    damping(T day_in, T factor_in)
    {

        day = day_in;
        factor = factor_in;

    }
};



/**
 * Paramters of the model(s):
 * T_inc (also sigma^(-1) or R_2^(-1)+R_3^(-1)): mean incubation period (default: 5.2);
 *          R_2^(-1) is the first part of the incubation time where the person is not yet infectioous
 *          R_3 is the exchange between asymptomatic carriers and infectious people; R_3^(-1) is the second part of the incubation time where the person is infectious WITHOUT showing symptoms
 * T_serint (also R_2^(-1)+0.5*R_3^(-1)): serial interval (default: 4.2);
 * T_infmild (also gamma^(-1) or R_4^(-1)): time a person remains infective after disease (if 'hospitalized' is considered a state, it does not apply to them but only to 'mildly infected' people in SECIR)
 * T_hosp2home (also R_5^(1)): duration for which the hospitalized patients not requiring further intensive care remain under general hospital care (=INF or R_5=0 in standard SEIR to waive influence of this parameter)
 * T_home2hosp (also R_6^(-1)): mean time a patient with mild symptoms spends at home before hospital admission due to worsening of the disease condition  (=INF or R_6=0 in standard SEIR to waive influence of this parameter)
 * T_hosp2icu (also R_7^(-1)): mean time a patient who entered the hospital will be hopistalized without ICU before being connected to an ICU  (=INF or R_7=0 in standard SEIR to waive influence of this parameter)
 * T_icu2home (also R_8^(-1)): mean time a patient is connected to an ICU before returning home (=INF or R_8=0 in standard SEIR to waive influence of this parameter)
 * T_infasy (also R_9^(-1)): mean time an asymptomatic person remains infective (=INF or R_9=0 in standard SEIR to waive influence of this parameter)
 * T_icu2death (also d; better would be R_10^(-1)): mean time a person needs ICU support before dying (=INF or R_10=0 in standard SEIR to waive influence of this parameter)
 * cont_freq (also R_1: contact frequency
 * alpha: share of asymptomatic cases
 * beta (Not the beta in SEIR model): risk of infection from the infected symptomatic patients
 * rho: H/I; hospitalized per infected (=0 in standard SEIR)
 * theta: U/H; intensive care units per hospitalized
 * delta: D/U; deaths per intensive care units
**/
template <typename T>
struct seirParam {

  int model;

  T base_reprod, b, cont_freq;
  T tinc_inv, tinfmild_inv;
  T tserint_inv, thosp2home_inv, thome2hosp_inv, thosp2icu_inv, ticu2home_inv, tinfasy_inv, ticu2death_inv;
  T alpha, beta, rho, theta, delta;

  // T nb_total, nb_exp, nb_car, nb_inf, nb_hosp, nb_icu, nb_rec, nb_dead;
  T nb_total_t0, nb_sus_t0, nb_exp_t0, nb_car_t0, nb_inf_t0, nb_hosp_t0, nb_icu_t0, nb_rec_t0, nb_dead_t0;

  // This defines a damping factor for a mitigation strategy for different points in time.
  std::vector<damping<T> > dampings;

  seirParam() {
    // assume an incubation period of 5.2 days;
    // an infectious period of (nonhospitalized) people (after disease) of 6 days
    // and a basis reproduction number (R0) of 2.7
    tinc_inv = 1.0/5.2; // 1.0/5.2 (in JS version)
    tinfmild_inv = 1.0/6.0; // 1.0/2.0 (in JS version)
    base_reprod = 2.7; // 3.5 (in JS version)

    // contact rate beta
    b = base_reprod*tinfmild_inv;

    model = 0; // standard is SEIR

    // Initial Population size
    nb_total_t0 = 10000;
    // Initial Number of exposed
    nb_exp_t0 = 50.0;
    // Intial Number of infected
    nb_inf_t0 = 10.0;
    // Initial Number of recovered
    nb_rec_t0 = 0.0;

    nb_sus_t0 = nb_total_t0-nb_exp_t0-nb_inf_t0-nb_rec_t0;

    // additional numbers for SECIR (SECIHURD) model
    nb_car_t0 = 0; nb_hosp_t0 = 0; nb_icu_t0 = 0; nb_dead_t0 = 0;
    // List of damping initially empty
    dampings.push_back(damping<T>(0.0, 1.0));
  }


  seirParam(T tinc, T tinfmild, T base_reprod_in, T nb_total_t0_in, T nb_exp_t0_in, T nb_inf_t0_in, T nb_rec_t0_in) 
  {
    tinc_inv = 1.0/tinc;
    tinfmild_inv = 1.0/tinfmild;
    base_reprod = base_reprod_in;

    // contact rate beta
    b = base_reprod*tinfmild_inv;

    model = 0; // standard is SEIR

    // Initial Population size
    nb_total_t0 = nb_total_t0_in;
    // Initial Number of exposed
    nb_exp_t0 = nb_exp_t0_in;
    // Intial Number of infected
    nb_inf_t0 = nb_inf_t0_in;
    // Initial Number of recovered
    nb_rec_t0 = nb_rec_t0_in;
    // List of damping initially empty
    dampings.push_back(damping<T>(0.0, 1.0));
  }

  seirParam(T tinc, T tinfmild, T tserint, T thosp2home, T thome2hosp, T thosp2icu, T ticu2home, T tinfasy, T ticu2death, T cont_freq_in, T alpha_in, T beta_in, T delta_in, T rho_in, T theta_in, 
                  T nb_total_t0_in, T nb_exp_t0_in, T nb_car_t0_in, T nb_inf_t0_in, T nb_hosp_t0_in, T nb_icu_t0_in, T nb_rec_t0_in, T nb_dead_t0_in) 
  {

    
    tinc_inv = 1.0/tinc; // inverse incubation period
    tinfmild_inv = 1.0/tinfmild; // inverse infectious period (of nonhospitalized cases)
    tserint_inv = 1.0/tserint; // inverse serial interval
    thosp2home_inv = 1.0/thosp2home; // inverse hospital to home time
    thome2hosp_inv = 1.0/thome2hosp; // inverse home to hospital time
    thosp2icu_inv = 1.0/thosp2icu; // inverse hospital to ICU time
    ticu2home_inv = 1.0/ticu2home; // inverse icu to home time
    tinfasy_inv = 1.0/tinfasy; // inverse infectious period for asymptomatic cases
    ticu2death_inv = 1.0/ticu2death; // inverse ICU to death time

    cont_freq = cont_freq_in; // contact frequency 
    alpha = alpha_in; // percentage of asymptomatic cases
    beta = beta_in; // risk of infection from the infected symptomatic patients
    rho = rho_in; // hospitalized per infected
    theta = theta_in; // icu per hospitalized
    delta = delta_in; // deaths per ICUs

    model = 1; 

    // Initial Population size
    nb_total_t0 = nb_total_t0_in;
    // Initial Number of exposed
    nb_exp_t0 = nb_exp_t0_in;
    // Intial Number of infected
    nb_inf_t0 = nb_inf_t0_in;
    // Initial Number of recovered
    nb_rec_t0 = nb_rec_t0_in;

    // Initial Number of carriers
    nb_car_t0 = nb_car_t0_in; 
    // Initial Number of hospitalized
    nb_hosp_t0 = nb_hosp_t0_in;
    // Initial Number of ICU
    nb_icu_t0 = nb_icu_t0_in; 
    // Initial Number of deaths
    nb_dead_t0 = nb_dead_t0_in;
    // List of damping initially empty
    dampings.push_back(damping<T>(0.0, 1.0));

    T dummy_R3 = 0.5/(tinfmild-tserint);

    nb_sus_t0 = nb_total_t0-nb_exp_t0-nb_car_t0-nb_inf_t0-nb_hosp_t0-nb_icu_t0-nb_rec_t0-nb_dead_t0;

    // only for output... (not explicitly used as single parameter in SECIR)
    base_reprod = cont_freq*((1-rho)*tinfmild_inv+dummy_R3*beta*(1-alpha)+rho*thome2hosp)/((dummy_R3*(1-alpha)+alpha*tinfasy_inv)*(tinfmild_inv*(1-rho)+rho*thome2hosp_inv))*nb_sus_t0/nb_total_t0;
  }

  void add_damping(const damping<T>& d) {
      dampings.push_back(d);
  }
    
};

/**
 * prints given parameters
 * @param[in] params the seirParam parameter object
 */
template <typename T>
void printSeirParams(struct seirParam<T> const &params)
{
    if(params.model == 0)
    {
        printf("\n SEIR model set.\n Parameters:\n\t Time incubation:\t %.4f \n\t Time infectious:\t %.4f \n\t b:\t %.4f \n\t N:\t %d \n\t E0:\t %d \n\t I0:\t %d \n\t R0:\t %d\n",
          1.0/params.tinc_inv, 1.0/params.tinfmild_inv, params.b, (int)params.nb_total_t0, (int)params.nb_exp_t0, (int)params.nb_inf_t0, (int)params.nb_rec_t0);
    }else{
        printf("\n SECIR (SECIHURD) model set.\n Parameters:\n\t Time incubation:\t %.4f \n\t Time infectious (mild):\t %.4f \n\t Serial interval:\t %.4f \n\t Time hosp.->home:\t %.4f \n\t Time home->hosp.:\t %.4f \n\t Time hosp.->icu:\t %.4f \n\t Time infectious (asymp.):\t %.4f \n\t Time icu->death:\t\t %.4f\n\t contact freq.:\t %.4f \n\t alpha:\t %.4f \n\t beta:\t %.4f \n\t delta:\t %.4f \n\t rho:\t %.4f \n\t theta:\t %.4f \n\t N0:\t %d \n\t E0:\t %d \n\t C0:\t %d\n\t I0:\t %d \n\t H0:\t %d \n\t U0:\t %d \n\t R0:\t %d \n\t D0:\t %d\n\t Calculated R_0: %.4f\n",
               1.0/params.tinc_inv, 1.0/params.tinfmild_inv, 1.0/params.tserint_inv,  1.0/params.thosp2home_inv, 1.0/params.thome2hosp_inv, 1.0/params.thosp2icu_inv, 1.0/params.tinfasy_inv, 1.0/params.ticu2death_inv,
            params.cont_freq, params.alpha, params.beta, params.delta, params.rho, params.theta,
            (int)params.nb_total_t0, (int)params.nb_exp_t0, (int)params.nb_car_t0, (int)params.nb_inf_t0, (int)params.nb_hosp_t0, (int)params.nb_icu_t0, (int)params.nb_rec_t0,  (int)params.nb_dead_t0,
               params.base_reprod);
    }
}


#endif // SEIR_PARAM
