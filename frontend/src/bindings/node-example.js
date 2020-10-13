var secirjs = require('./secirjs.js');

function toArray(vector) {
  var arr = [];
  for (var i = 0; i < vector.size(); ++i) {
    arr.push(vector.get(i));
  }
  return arr;
}

secirjs().then(function (secir) {
  var times = new secir.StageTimes();
  times.set_incubation(5.2); // R_2 ^ (-1) + R_3 ^ (-1)
  times.set_infectious_mild(6); // 4 - 14(=R4 ^ (-1))
  times.set_serialinterval(4.2); // 4 - 4.4 // R_2^(-1)+0.5*R_3^(-1)
  times.set_hospitalized_to_home(12); // 7 - 16(=R5 ^ (-1))
  times.set_home_to_hospitalized(5); // 2.5 - 7(=R6 ^ (-1))
  times.set_hospitalized_to_icu(2); // 1 - 3.5(=R7 ^ (-1))
  times.set_icu_to_home(8); // 5 - 16(=R8 ^ (-1))
  times.set_infectious_asymp(6.2); // (= R9 ^ (-1)= R_3 ^ (-1) + 0.5 * R_4 ^ (-1))
  times.set_icu_to_death(5); // 3.5 - 7(=R5 ^ (-1))

  var probs = new secir.Probabilities();
  probs.set_infection_from_contact(0.05);
  probs.set_carrier_infectability(0.67);
  probs.set_asymp_per_infectious(0.09); // 0.01 - 0.16
  probs.set_risk_from_symptomatic(0.25); // 0.05 - 0.5
  probs.set_hospitalized_per_infectious(0.2); // 0.1 - 0.35
  probs.set_icu_per_hospitalized(0.25); // 0.15 - 0.4
  probs.set_dead_per_icu(0.3); // 0.15 - 0.77

  var people = new secir.Populations();
  people.set_total_t0(10000);
  people.set_exposed_t0(100);
  people.set_carrier_t0(50);
  people.set_infectious_t0(50);
  people.set_hospital_t0(20);
  people.set_icu_t0(10);
  people.set_recovered_t0(10);
  people.set_dead_t0(0);

  // set the params required fseciror the simulation
  var param = new secir.SecirParams();
  param.times = times;
  param.probabilities = probs;
  param.populations = people;
  times.delete();
  probs.delete();
  people.delete();

  var params = new secir.VectorSecirParams();
  params.push_back(param);
  param.delete();

  secir.print_secir_params(params);

  params[0].get_cont_freq_matrix().set_cont_freq(10, 0, 0); // 0.2 - 0.75

  // emulate some mitigations
  var d1 = new secir.Damping(23, 0.8);
  var d2 = new secir.Damping(25, 0.75);
  var d3 = new secir.Damping(27, 0.7);
  params[0].get_cont_freq_matrix().add_damping(d1, 0, 0);
  params[0].get_cont_freq_matrix().add_damping(d2, 0, 0);
  params[0].get_cont_freq_matrix().add_damping(d3, 0, 0);
  d1.delete();
  d2.delete();
  d3.delete();

  result = secir.simulate(0, 100, 0.1, params);
  params.delete();

  console.log(toArray(result.t));

  result.delete();
});
