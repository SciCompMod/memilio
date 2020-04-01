export {
    makeSeirParam,
    SeirParam,
    simulate_seir
};

import { euler } from './euler.js';

class SeirParam {
    constructor() {
      // Assume an incubation period of 5.2 days
      this.a = 0.0192;
      // contact rate beta
      this.b = 1.75;
      // Assume infectious period of 2 days
      this.g = 0.5;
      // Initial Population size
      this.N = 10000;
      // Initial Number of exposed
      this.E0 = 1000.;
    }
  }
  
  // SeirParam Factory
  function makeSeirParam()
  {
    return new SeirParam();
  }
  
  /**
   * Computes the current time-derivative of the seir model
   * @param {} y Current SEIR values at t
   * @param {*} params SEIR Model parameters, created by seir_param
   */
  function seir(y, params) {
    const [S,E,I,R] = y;
    
    var dydt = [
      -params.b*S*I/params.N,
       params.b*S*I/params.N - params.a*E,
       params.a*E - params.g*I,
       params.g*I
  
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
  function simulate_seir(t0, tmax, dt, params)
  {
    const seir_0 = [params.N-params.E0, params.E0, 0.0, 0.0];
  
    const n = Math.ceil((tmax-t0) / dt);
    const t = tf.linspace(t0, tmax, n);
  
    // Euler integration - might switch to RK4 though
    var res = euler(function(t, y) {
      return seir(y, params);
    }, seir_0, t);
    const [S, E, I, R] = tf.split(res, 4, 0);
  
    return {
      t: t.arraySync(),
      S: S.arraySync()[0],
      E: E.arraySync()[0],
      I: I.arraySync()[0],
      R: R.arraySync()[0],
    };
  }
  