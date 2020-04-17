#ifndef SECIR_H
#define SECIR_H

#include <vector>
#include <cmath>
#include <cassert>
#include <cstdio>

// #include <epidemiology/euler.h>
#include <epidemiology/adapt_rk.h>
// #include <epidemiology/seirParam.h>


/**
 * Returns the damping factor
 *
 * @param[in] damping_array Array of dampings
 * @param[in] t Current day
 */
template <typename T>
T getDampingFactor(std::vector<struct damping<T> > const& damping_array, T day)
{
    // we assume, that the data_array is ordered in ascending order
    size_t ilow = 0;
    size_t ihigh = damping_array.size() - 1;

    // check extrapolation cases
    if (day < damping_array[ilow].day) {
        return damping_array[ilow].factor;
    }

    if (day >= damping_array[ihigh].day) {
        return damping_array[ihigh].factor;
    }

    // now do the search
    while (ilow < ihigh - 1) {
        size_t imid = (ilow + ihigh) / 2;
        if (damping_array[ilow].day <= day && day < damping_array[imid].day) {
            ihigh = imid;
        }
        else if(damping_array[imid].day <= day && day < damping_array[ihigh].day) {
            ilow = imid;
        }
        else {
            // this case can only occur, if
            // input data are not ordered
            return 1e16;
        }
    }

    assert(damping_array[ilow].day <= day && day < damping_array[ilow + 1].day);
    return damping_array[ilow].factor;
}

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
  

/**
 * Computes the current time-derivative of S, E, I, and R in the SEIR model
 * @tparam T the datatype of the cases
 * @param[in] params SEIR Model parameters, created by seir_param
 * @param[in] y current  S, E, I, and R values at t; y: [0:S, 1:E, 2:I, 3:R]
 * @param[in] t time / current day
 * @param[out] dydt the values of the time derivatices of S, E, I, and R
 */
template <typename T>
void seir_getDerivatives(struct seirParam<T> const &params, std::vector<T> const &y, const T t, std::vector<T> &dydt)
{

  T b_eff = params.b * getDampingFactor(params.dampings, t);
  T divN = 1.0/params.nb_total_t0;

  dydt[0] = -b_eff * y[0] * y[2] * divN;
  dydt[1] =  b_eff * y[0] * y[2] * divN - params.tinc_inv * y[1];
  dydt[2] =  params.tinc_inv * y[1] - params.tinfmild_inv * y[2];
  dydt[3] =  params.tinfmild_inv * y[2];

}


/**
 * Computes the current time-derivative of S, E, I, and R in the SEIR model
 * @param[in] params SEIR Model parameters, created by seir_param
 * @tparam T the datatype of the cases
 * @param[in] y current  S, E, I, and R values at t; y: [0:S, 1:E, 2:I, 3:R]
 * @param[in] t time / current day
 * @param[out] dydt the values of the time derivatices of S, E, I, and R
 */
template <typename T>
void secir_getDerivatives(struct seirParam<T> const &params, std::vector<T> const &y, const T t, std::vector<T> &dydt)
{
  
  // 0: S,      1: E,     2: C,     3: I,     4: H,     5: U,     6: R,     7: D
  T cont_freq_eff = params.cont_freq * getDampingFactor(params.dampings, t);
  T divN = 1.0/params.nb_total_t0;

  T dummy_S = cont_freq_eff * y[0] * divN * (y[2]+params.beta*y[3]);

  T dummy_R2 = 1.0/(2*(1.0/params.tserint_inv)-(1.0/params.tinfmild_inv)); // R2 = 1/(2SI-IP)
  T dummy_R3 = 0.5/((1.0/params.tinfmild_inv)-(1.0/params.tserint_inv)); // R3 = 1/(2(IP-SI))

  dydt[0] = -dummy_S; // -R1*(C+beta*I)*S/N0
  dydt[1] =  dummy_S - dummy_R2*y[1]; // R1*(C+beta*I)*S/N0-R2*E - R2*E
  dydt[2] =  dummy_R2*y[1] - ((1-params.alpha)*dummy_R3+params.alpha*params.tinfasy_inv)*y[2];
  dydt[3] =  (1-params.alpha)*dummy_R3*y[2] - ((1-params.rho)*params.tinfmild_inv+params.rho*params.thome2hosp_inv)*y[3];
  dydt[4] = params.rho*params.thome2hosp_inv*y[3] - ((1-params.theta)*params.thosp2home_inv + params.theta*params.thosp2icu_inv)*y[4];
  dydt[5] = params.theta*params.thosp2icu_inv*y[4] - ((1-params.delta)*params.ticu2home_inv+params.delta*params.ticu2death_inv)*y[5];
  dydt[6] = params.alpha*params.tinfasy_inv*y[2] + (1-params.rho)*params.tinfmild_inv + (1-params.theta)*params.thosp2home_inv*y[4] + (1-params.delta)*params.ticu2home_inv*y[5];
  dydt[7] = params.delta*params.ticu2death_inv*y[5];

}
  
/**
 * Computes the values of S, E, I, and R at the current time step in the SEIR model
 * by using a numerical integration scheme (only explicit Euler at the moment)
 * @tparam T the datatype of the cases
 * @param[in] dydt the values of the time derivatives of S, E, I, and R
 * @param[in] dt the current time step
 * @param[in] y the vector of  S, E, I, and R values at t
 * @param[out] result the result vector of  S, E, I, and R values at t+1
 */
template <typename T>
void integrate_euler(const std::vector<T>& y, std::vector<T> const& dydt, const T dt, std::vector<T>& result)
{
    explicit_euler(y, dydt, dt, result);
}


/**
 * Computes the values of S, E, I, and R at the current time step in the SEIR model
 * by using an adaptive time stepping numerical integration scheme 
 * @tparam T the datatype of the cases
 * @param[in] tab the upper part of the Butcher teableau for the adaptive scheme
 * @param[in] tab_final the final row(s) (two rows) for the two different integration schemes of different order
 * @param[in] params SEIR Model parameters, created by seir_param
 * @param[in] y the vector of  S, E, I, and R values at t
 * @param[in] f the values of the time derivatives of S, E, I, and R
 * @param[in] tol_abs the required absolute tolerance for the comparison with the Fehlberg approximation
 * @param[in] tol_rel the required relative tolerance for the comparison with the Fehlberg approximation
 * @param[in] dtmin the minimum step size
 * @param[in] dtmax the maximum step size
 * @param[inout] t current time step h=dt
 * @param[inout] dt current time step h=dt
 * @param[out] result the result vector of  S, E, I, and R values at t+1
 */
template <typename T, typename rhsDerivatives>
bool integrate_ark(struct tableau<T> const &tab, struct tableau_final<T> const &tab_final, struct seirParam<T> const& params, const std::vector<T>& y, rhsDerivatives const &f, const T tol_abs, const T tol_rel, const T dtmin, const T dtmax, T &t, T &dt, std::vector<T>& result)
{
  bool failed_step_size_adapt = false;
    
  failed_step_size_adapt = adapt_rk(tab, tab_final, params, y, f, tol_abs, tol_rel, dtmin, dtmax, t, dt, result);

  return failed_step_size_adapt;
}


/**
 * Computes the seir curve by integration
 * @param[in] seir_0 Initial S, E, I, and R values at t0
 * @param[in] t0 start time of simulation
 * @param[in] tmax end time of simulation
 * @param[in] dt initial time step
 * @param[in] params SEIR model parameters
 */
template <typename T>
void simulate(const double t0, const double tmax, T dt, struct seirParam<T> const &params, std::vector<std::vector<T> > &seir)
{
  size_t nb_steps = (int)(ceil((tmax-t0) / dt)); // estimated number of time steps (if equidistant)

  std::vector<T> dydt;  

  if(params.model == 0)
  {
    seir = std::vector<std::vector<T>>(nb_steps+1, std::vector<T>(4,(T)0)); // prepare memory for equidistant step size
    
    //initial conditions
    seir[0][0] = params.nb_sus_t0;  
    seir[0][1] = params.nb_exp_t0; 
    seir[0][2] = params.nb_inf_t0; 
    seir[0][3] = params.nb_rec_t0;
    
    dydt = std::vector<T>(4,0);

  }else{
    seir = std::vector<std::vector<T>>(nb_steps+1, std::vector<T>(8,(T)0)); // prepare memory for equidistant step size
    
    //initial conditions
    seir[0][0] = params.nb_sus_t0;
    seir[0][1] = params.nb_exp_t0; 
    seir[0][2] = params.nb_car_t0; 
    seir[0][3] = params.nb_inf_t0;
    seir[0][4] = params.nb_hosp_t0;
    seir[0][5] = params.nb_icu_t0;
    seir[0][6] = params.nb_rec_t0;
    seir[0][7] = params.nb_dead_t0;
    
    dydt = std::vector<T>(8,0);
  }

  #ifdef ARK_H
  tableau<double> rkf45_tableau = tableau<double>();
  tableau_final<double> rkf45_tableau_final = tableau_final<double>();
  double tol_abs = 1e-1;
  double tol_rel = 1e-4;
  double dtmin = 1e-3;
  double dtmax = 1;  
  bool failed_step_size_adapt = false;
  #endif

  T t = t0;
  T t_prev = t0;
  size_t i = 0;
  while(t_prev < tmax)
  {
    if(t > tmax) // possible for adaptive step size
    {
      dt = tmax - t_prev;
      if(dt < 0.1*dtmin)
      {
        break;
      }
    }
    t_prev = t;
    if(t > tmax) // possible for adaptive step size
    {
      t = tmax;
    }

    // printf("%d\t", (int)seir[i][0]);
    #ifdef EULER_H
    if(params.model == 0)
    {
      seir_getDerivatives(params, seir[i], t, dydt);
    }else{
      secir_getDerivatives(params, seir[i], t, dydt);
    }
    #endif

    #ifdef ARK_H
      bool failed_step_size_adapt_dummy = false;
      if(params.model == 0)
      {
        if(i+1 >= seir.size())
        {
            std::vector<std::vector<T>> vecAppend(20, std::vector<T>(4,(T)0));
            seir.insert(seir.end(), vecAppend.begin(), vecAppend.end());
        }
        failed_step_size_adapt_dummy = integrate_ark(rkf45_tableau, rkf45_tableau_final, params, seir[i], seir_getDerivatives<T>, tol_abs, tol_rel, dtmin, dtmax, t, dt, seir[i+1]);
      }else{
        if(i+1 >= seir.size())
        {
            std::vector<std::vector<T>> vecAppend(20, std::vector<T>(8,(T)0));
            seir.insert(seir.end(), vecAppend.begin(), vecAppend.end());
        }
        failed_step_size_adapt_dummy = integrate_ark(rkf45_tableau, rkf45_tableau_final, params, seir[i], secir_getDerivatives<T>, tol_abs, tol_rel, dtmin, dtmax, t, dt, seir[i+1]);
      }
      if(failed_step_size_adapt_dummy)
      {
        failed_step_size_adapt = true;
      }
    #elif EULER_H
      integrate_euler(seir[i], dydt, dt, seir[i+1]);
      t += dt;
    #endif
    
    i++;
  }

  // cut empty elements (makes more sense for adaptive time step size)
  if(seir.size() > i)
  {
    seir.resize(i);
  }

  if(!failed_step_size_adapt)
  {
    printf("\n Adaptive step sizing successful to tolerances.\n");
  }else{
    printf("\n Adaptive step sizing failed.");
  }

}
#endif // SECIR_H
