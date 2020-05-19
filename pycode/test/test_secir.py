import unittest
from epidemiology.secir import ContactFrequencyMatrix, Damping, SecirParams, simulate, StageTimes, Probabilities, Populations

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
        times.set_infectious_asymp(6.2)  # (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
        times.set_icu_to_death(5.)  # 3.5-7 (=R5^(-1))

        probs = Probabilities()
        probs.set_asymp_per_infectious(0.09)  # 0.01-0.16
        probs.set_risk_from_symptomatic(0.25)  # 0.05-0.5
        probs.set_hospitalized_per_infectious(0.2)  # 0.1-0.35
        probs.set_icu_per_hospitalized(0.25)  # 0.15-0.4
        probs.set_dead_per_icu(0.3)  # 0.15-0.77

        people = Populations()
        people.set_total_t0(10000)
        people.set_exposed_t0(100)
        people.set_carrier_t0(50)
        people.set_infectious_t0(50)
        people.set_hospital_t0(20)
        people.set_icu_t0(10)
        people.set_recovered_t0(10)
        people.set_dead_t0(0)

        # set the params required for the simulation
        params = [SecirParams()]
        params[0].times = times
        params[0].probabilities = probs
        params[0].populations = people

        cont_freq_matrix = ContactFrequencyMatrix()
        cont_freq_matrix.set_cont_freq(0.5, 0, 0)  # 0.2-0.75

        self.params = params
        self.cont_freq_matrix = cont_freq_matrix    

    def test_simulate_simple(self):
      result = simulate(t0=0., tmax=100., dt=0.1, cont_freq_matrix=self.cont_freq_matrix, params=self.params)
      self.assertAlmostEqual(result.t[0], 0.)
      self.assertAlmostEqual(result.t[1], 0.1)
      self.assertAlmostEqual(result.t[-1], 100.)

if __name__ == '__main__':
    unittest.main()