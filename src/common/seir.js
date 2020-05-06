import * as tf from '@tensorflow/tfjs';
import { euler } from './euler.js';

/**
 * This defined a damping factor for a
 * mitigation strategy for one point in time.
 */
class Damping {
  constructor(day, factor) {
    this.day = day;
    this.factor = factor;
  }
}

/**
 *  Finds index idx such that
 *
 *    get_value_function(data_array[idx]) <= t < get_value_function(data_array[idx+1])
 *
 *  We are using a bracketing approach
 */
function bracket(data_array, value, get_value_function) {
  // we assume, that the data_array is ordered in ascending order
  let ilow = 0;
  let ihigh = data_array.length - 1;

  // check extrapolation cases
  if (value < get_value_function(data_array[ilow])) {
    return ilow;
  }

  if (value >= get_value_function(data_array[ihigh])) {
    return ihigh;
  }

  // now do the search
  while (ilow < ihigh - 1) {
    let imid = Math.floor((ilow + ihigh) / 2);
    if (
      get_value_function(data_array[ilow]) <= value &&
      value < get_value_function(data_array[imid])
    ) {
      ihigh = imid;
    } else if (
      get_value_function(data_array[imid]) <= value &&
      value < get_value_function(data_array[ihigh])
    ) {
      ilow = imid;
    } else {
      // this case can only occur, if
      // input data are not ordered
      return data_array.length;
    }
  }

  return ilow;
}

/**
 * Returns the damping factor rho(t)
 *
 * @param {*} damping_array Array of dampings
 * @param {*} t Current day
 */
function getDampingFactor(damping_array, t) {
  if (damping_array.length === 0) {
    return 1;
  } else if (damping_array.length === 1) {
    return t < damping_array[0].day ? 1 : damping_array[0].factor;
  }

  // find index i such that  dampings[idx].day <= t < dampings[idx+1].day
  let idx = bracket(damping_array, t, function (damping) {
    return damping.day;
  });

  return damping_array[idx].factor;
}

class SeirParam {
  constructor() {
    // Assume an incubation period of 5.2 days, a = 1/t_incubation
    this.a = 0.192;
    // contact rate beta
    this.b = 1.75;
    // Assume infectious period of 2 days
    this.g = 0.5;
    // Initial Population size
    this.N = 10000;
    // Initial Number of exposed
    this.E0 = 10000.0;
    // Intial Number of infected
    this.I0 = 10000.0;
    // Initial Number of recovered
    this.R0 = 10000.0;
    // List of "Damping" objects / rhos that can be used
    // to model social distancing
    this.dampings = [];
  }
}

// SeirParam Factory
function makeSeirParam() {
  return new SeirParam();
}

/**
 * Computes the current time-derivative of the seir model
 * @param {} y Current SEIR values at t
 * @param {} t Time / Current day
 * @param {*} params SEIR Model parameters, created by seir_param
 */
function seir(y, t, params) {
  const [S, E, I] = y;

  var b_eff = params.b * getDampingFactor(params.dampings, t);

  var dydt = [
    (-b_eff * S * I) / params.N,
    (b_eff * S * I) / params.N - params.a * E,
    params.a * E - params.g * I,
    params.g * I
  ];
  return dydt;
}

/**
 * Computes the seir curve by integration
 * @param {} seir_0 Initial S/E/I/R Numbers at t0
 * @param {*} t0 Start time
 * @param {*} tmax  End time of simulation
 * @param {*} dt Time step
 * @param {*} params SEIR model parameters
 */
function simulate_seir(t0, tmax, dt, params) {
  //initial conditions
  const seir_0 = [
    params.N - params.E0 - params.I0 - params.R0,
    params.E0,
    params.I0,
    params.R0
  ];

  const n = Math.ceil((tmax - t0) / dt);
  const t = tf.linspace(t0, tmax, n);

  // Euler integration - might switch to RK4 though
  var res = euler(
    function (t, y) {
      return seir(y, t, params);
    },
    seir_0,
    t
  );
  const [S, E, I, R] = tf.split(res, 4, 0);

  return {
    t: t.arraySync(),
    S: S.arraySync()[0],
    E: E.arraySync()[0],
    I: I.arraySync()[0],
    R: R.arraySync()[0]
  };
}

export { makeSeirParam, SeirParam, simulate_seir, Damping };
