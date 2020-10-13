import unittest
from epidemiology.secir import (UncertainContactMatrix, ContactFrequencyMatrix, Damping, SecirParams,
                                simulate, StageTimes, Probabilities, Populations, SecirCompartments, SecirSimulation)

class Test_secir_integration(unittest.TestCase):

    def setUp(self):
        times = StageTimes()
        times.set_incubation(5.2)  # R_2^(-1)+R_3^(-1)
        times.set_infectious_mild(6.)  # 4-14  (=R4^(-1))
        times.set_serialinterval(4.2)   # 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        times.set_hospitalized_to_home(12.)  # 7-16 (=R5^(-1))
        times.set_home_to_hospitalized(5.)  # 2.5-7 (=R6^(-1))
        times.set_hospitalized_to_icu(2.)  # 1-3.5 (=R7^(-1))
        times.set_icu_to_home(8.)  # 5-16 (=R8^(-1))
        times.set_icu_to_death(5.)  # 3.5-7 (=R5^(-1))

        probs = Probabilities()
        probs.set_infection_from_contact(1.0)
        probs.set_asymp_per_infectious(0.09)  # 0.01-0.16
        probs.set_risk_from_symptomatic(0.25)  # 0.05-0.5
        probs.set_hospitalized_per_infectious(0.2)  # 0.1-0.35
        probs.set_icu_per_hospitalized(0.25)  # 0.15-0.4
        probs.set_dead_per_icu(0.3)  # 0.15-0.77

        people = Populations([1, SecirCompartments.SecirCount])
        people.set([0, SecirCompartments.S], 7600)
        people.set([0, SecirCompartments.E], 100)
        people.set([0, SecirCompartments.C], 50)
        people.set([0, SecirCompartments.I], 50)
        people.set([0, SecirCompartments.H], 20)
        people.set([0, SecirCompartments.U], 10)
        people.set([0, SecirCompartments.R], 10)
        people.set([0, SecirCompartments.D], 0)

        # set the params required for the simulation
        params = SecirParams(1)
        params.times = [times]
        params.probabilities = [probs]
        params.populations = people

        params.get_contact_patterns().get_cont_freq_mat().set_cont_freq(0.5, 0, 0)

        params.apply_constraints()

        self.params = params

    def test_simulate_simple(self):
      result = simulate(t0=0., tmax=100., dt=0.1, params=self.params)
      self.assertAlmostEqual(result.get_time(0), 0.)
      self.assertAlmostEqual(result.get_time(1), 0.1)
      self.assertAlmostEqual(result.get_last_time(), 100.)

    def test_simulation_simple(self):
        sim = SecirSimulation(self.params, t0 = 0., dt = 0.1)
        sim.advance(tmax = 100.)
        result = sim.result
        self.assertAlmostEqual(result.get_time(0), 0.)
        self.assertAlmostEqual(result.get_time(1), 0.1)
        self.assertAlmostEqual(result.get_last_time(), 100.)

if __name__ == '__main__':
    unittest.main()
