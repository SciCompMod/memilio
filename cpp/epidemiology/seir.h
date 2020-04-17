#ifndef SEIR_H
#define SEIR_H

#include <vector>
#include <cmath>
#include <cassert>

#include <epidemiology/euler.h>


/**
 * This defined a damping factor for a
 * mitigation strategy for one point in time.
 */
template <typename T>
struct damping {
    T day;
    T factor;

    damping(int day_in, T factor_in)
    {

        day = day_in;
        factor = factor_in;

    }
};


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


template <typename T>
struct seirParam {

    T a, b, g, N, E0, I0, R0;

    // This defines a damping factor for a mitigation strategy for different points in time.
    std::vector<damping<T> > dampings;

    seirParam()
    {
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
        dampings.push_back(damping<T>(0, 1.0));
    }


    seirParam(T a_in, T b_in, T g_in, T N_in, T E0_in, T I0_in, T R0_in)
    {
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
        // List of damping initially
        dampings.push_back(damping<T>(0, 1.0));
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
void printSeirParams(struct seirParam<T> const& params)
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
void seir_getDerivatives(std::vector<T> const& y, struct seirParam<T> const& params, const T t, std::vector<T>& dydt)
{

    T b_eff = params.b * getDampingFactor(params.dampings, t);
    T divN = 1.0 / params.N;

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
void seir_integrate(std::vector<T> const& dydt, struct seirParam<T> const& params, const T dt, const std::vector<T>& y, std::vector<T>& result)
{
    explicit_euler(y, dt, dydt, result);
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
void simulate_seir(const T t0, const T tmax, const T dt, struct seirParam<T> const& params, std::vector<std::vector<T> >& seir)
{
    size_t nb_steps = (int)(ceil((tmax - t0) / dt)); // estimated number of time steps (if equidistant)

    seir = std::vector<std::vector<T> >(nb_steps, std::vector<T>(4, (T)0)); // prepare memory for equidistant step size

    //initial conditions
    seir[0][0] = params.N - params.E0 - params.I0 - params.R0;
    seir[0][1] = params.E0;
    seir[0][2] = params.I0;
    seir[0][3] = params.R0;

    std::vector<T> dydt(4, 0);

    T t = t0;
    for (size_t i = 0; i < nb_steps - 1; ++i) {

        t = t0 + i*dt;

        seir_getDerivatives(seir[i], params, t, dydt);

        seir_integrate(dydt, params, dt, seir[i], seir[i+1]);
    }
}

#endif // SEIR_H
