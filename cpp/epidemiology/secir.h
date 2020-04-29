#ifndef SECIR_H
#define SECIR_H

#include <vector>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <algorithm>

//#include <epidemiology/euler.h>
#include <epidemiology/adapt_rk.h>
#include <epidemiology/seirParam.h>


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
 * Computes the seir curve by integration
 * @param[in] seir_0 Initial S, E, I, and R values at t0
 * @param[in] t0 start time of simulation
 * @param[in] tmax end time of simulation
 * @param[in] dt initial time step
 * @param[in] params SEIR model parameters
 *
 * @returns Vector of times t
 */
template <typename T>
std::vector<T> simulate(const double t0, const double tmax, T dt, struct seirParam<T> const &params, std::vector<std::vector<T> > &seir)
{
  size_t nb_steps = (int)(ceil((tmax-t0) / dt)); // estimated number of time steps (if equidistant)

  std::vector<T> dydt;

  size_t n_params = params.model == 0 ? 4 : 8;
  seir = std::vector<std::vector<T>>(nb_steps+1, std::vector<T>(n_params,(T)0)); // prepare memory for equidistant step size
  std::vector<T> vec_times(nb_steps+1, 0.);

  dydt = std::vector<T>(n_params,0);

  if(params.model == 0)
  {
    //initial conditions
    seir[0][0] = params.nb_sus_t0;  
    seir[0][1] = params.nb_exp_t0; 
    seir[0][2] = params.nb_inf_t0; 
    seir[0][3] = params.nb_rec_t0;
  }
  else{
    //initial conditions
    seir[0][0] = params.nb_sus_t0;
    seir[0][1] = params.nb_exp_t0; 
    seir[0][2] = params.nb_car_t0; 
    seir[0][3] = params.nb_inf_t0;
    seir[0][4] = params.nb_hosp_t0;
    seir[0][5] = params.nb_icu_t0;
    seir[0][6] = params.nb_rec_t0;
    seir[0][7] = params.nb_dead_t0;
  }

  auto secir_fun = [&params](std::vector<T> const &y, const T t, std::vector<T> &dydt) {
      return secir_getDerivatives<T>(params, y, t, dydt);
  };

#ifdef ARK_H
  T dtmin = 1e-3;
  T dtmax = 1.;
  RKIntegrator<T> integrator(secir_fun, dtmin, dtmax);
  integrator.set_abs_tolerance(1e-1);
  integrator.set_rel_tolerance(1e-4);
#endif
#ifdef EULER_H
  T dtmin = dt;
  T dtmax = dt;
  EulerIntegrator<T> integrator(secir_fun);
#endif

  bool step_okay = true;

  T t = t0;
  T t_prev = t0;
  size_t i = 0;
  vec_times[0] = t0;
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
    t = std::min(t, tmax); // possible for adaptive step size

    // printf("%d\t", (int)seir[i][0]);

    if(i+1 >= seir.size())
    {
        std::vector<std::vector<T>> vecAppend(20, std::vector<T>(n_params,(T)0));
        seir.insert(seir.end(), vecAppend.begin(), vecAppend.end());
        vec_times.resize(vec_times.size() + 20, 0.);
    }
    step_okay = integrator.step(seir[i], t, dt, seir[i+1]);
    vec_times[i+1] = t;

    i++;
  }

  // cut empty elements (makes more sense for adaptive time step size)
  if(seir.size() > i)
  {
    seir.resize(i);
    vec_times.resize(i);
  }

  if(step_okay)
  {
    printf("\n Adaptive step sizing successful to tolerances.\n");
  }else{
    printf("\n Adaptive step sizing failed.");
  }

  return vec_times;

}
#endif // SECIR_H
