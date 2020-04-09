#include <stdlib.h>
#include <vector>
#include <unordered_map>
#include <math.h> 
#include "euler.h"

/**
 * Returns the damping factor rho(t)
 * 
 * @tparam T the datatype of the cases
 * @param[in] damping_map Hash map of dampings
 * @param[in] t Current day
 */
template <typename T>
T getDampingFactor(std::unordered_map<int, T> const &damping_map, int const &t) {

  typename std::unordered_map<int, T>::const_iterator map_iter = damping_map.find(t); 

  if(map_iter == damping_map.end())
  {
    return 1;
  }else{
    return (*map_iter).second;
  }

}


template <typename T>
struct seirParam {

  T a,b,g,N,E0,I0,R0;

  // This defines a damping factor for a mitigation strategy for different points in time.
  std::unordered_map<int, T> dampings;

  seirParam() {
    // Assume an incubation period of 5.2 days, a = 1/t_incubation
    a = 0.192;
    // contact rate beta
    b = 1.75;
    // Assume infectious period of 2 days
    g = 0.5;
    // Initial Population size
    N = 10000;
    // Initial Number of exposed
    E0 = 50.0;
    // Intial Number of infected
    I0 = 10.0;
    // Initial Number of recovered
    R0 = 0.0;
    // List of damping initially empty
    // ...
  }


    seirParam(T a_in, T b_in, T g_in, T N_in, T E0_in, T I0_in, T R0_in) {
    // Assume an incubation period of t_incubation = 1/a_in
    a = a_in;
    // contact rate beta
    b = b_in;
    // Assume infectious period of 1/g_in days
    g = g_in;
    // Initial Population size
    N = N_in;
    // Initial Number of exposed
    E0 = E0_in;
    // Intial Number of infected
    I0 = I0_in;
    // Initial Number of recovered
    R0 = R0_in;
    // List of damping initially empty
    // ...
  }

  ~seirParam()
  {
    dampings.clear();
  }
    
};

/**
 * prints given parameters
 * @param[in] params the seirParam parameter object
 */
template <typename T>
void printSeirParams(struct seirParam<T> const &params)
{
  printf("\n Parameters set:\n\t a:\t %.4f \n\t b:\t %.4f \n\t g:\t %.4f \n\t N:\t %d \n\t E0:\t %d \n\t I0:\t %d \n\t R0:\t %d\n", params.a, params.b, params.g, (int)params.N, (int)params.E0, (int)params.I0, (int)params.R0);
}
  

/**
 * Computes the current time-derivative of S, E, I, and R in the SEIR model
 * @tparam T the datatype of the cases
 * @param[in] y current  S, E, I, and R values at t; y: [0:S, 1:E, 2:I, 3:R]
 * @param[in] params SEIR Model parameters, created by seir_param
 * @param[in] t time / current day
 * @param[out] dydt the values of the time derivatices of S, E, I, and R
 */
template <typename T>
void seir_getDerivatives(std::vector<T> const &y, struct seirParam<T> const &params, const int t, std::vector<T> &dydt)
{

  T b_eff = params.b * getDampingFactor(params.dampings, t);
  T divN = 1.0/params.N;

  dydt[0] = -b_eff * y[0] * y[2] * divN;
  dydt[1] =  b_eff * y[0] * y[2] * divN - params.a * y[1];
  dydt[2] =  params.a * y[1] - params.g * y[2];
  dydt[3] =  params.g * y[2];

}
  
/**
 * Computes the values of S, E, I, and R at the current time step in the SEIR model
 * by using a numerical integration scheme (only explicit Euler at the moment)
 * @tparam T the datatype of the cases
 * @param[out] dydt the values of the time derivatices of S, E, I, and R
 * @param[in] dt the current time step
 * @param[in] params SEIR Model parameters, created by seir_param
 * @param[in] i the position of y(t) in y[0...n]
 * @param[inout] y the vector of  S, E, I, and R values from 0 to t, yet empty for t+1 to n
 */
template <typename T>
void seir_integrate(std::vector<T> const &dydt, struct seirParam<T> const &params, const T dt, const size_t i, std::vector<std::vector<T> > &y)
{

  explicit_euler(y[i], dt, dydt, y[i+1]);

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
void simulate_seir(const double t0, const double tmax, const T dt, struct seirParam<T> const &params, std::vector<std::vector<T> > &seir)
{
  size_t nb_steps = (int)(ceil((tmax-t0) / dt)); // estimated number of time steps (if equidistant)

  seir = std::vector<std::vector<T> >(nb_steps+1, std::vector<T>(4,(T)0)); // prepare memory for equidistant step size

  //initial conditions
  seir[0][0] = params.N-params.E0-params.I0-params.R0;  
  seir[0][1] = params.E0; 
  seir[0][2] = params.I0; 
  seir[0][3] = params.R0;

  std::vector<T> dydt(4,0);

  T t = t0;
  size_t i = 0;
  while(tmax-t > 1e-7)
  {

    printf("%d\t", (int)seir[i][0]);

    seir_getDerivatives(seir[i], params, t, dydt);

    seir_integrate(dydt, params, dt, i, seir);

    t += dt;
    i++;
  }

  // cut empty elements (makes more sense for adaptive time step size)
  if(seir.size() > i)
  {
    seir.resize(i);
  }


}

int main()
{
  double t0 = 0;
  double tmax = 200;
  double dt = 0.1;

  seirParam<double> params = seirParam<double>();

  printSeirParams(params);

  std::vector<std::vector<double> > seir(0);

  simulate_seir(t0, tmax, dt, params, seir);

}
