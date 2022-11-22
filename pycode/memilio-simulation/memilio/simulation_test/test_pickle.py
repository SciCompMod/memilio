import pickle
import unittest

import numpy as np

import memilio.simulation as msim


class Test_Pickle(unittest.TestCase):
    def test_date(self):
        test = msim.Date(1, 2, 3)

        data = pickle.dumps(test)
        pickle_test = pickle.loads(data)

        self.assertEqual(pickle_test.year, 1)
        self.assertEqual(pickle_test.month, 2)
        self.assertEqual(pickle_test.day, 3)

    def test_uncertain_value(self):
        test = msim.UncertainValue(2.2)

        data = pickle.dumps(test)
        pickle_test = pickle.loads(data)

        self.assertEqual(pickle_test.value, 2.2)

    def test_distribution(self):

        test = msim.ParameterDistributionNormal(0, 1, 0.4, 0.1)

        data = pickle.dumps(test)
        pickle_test = pickle.loads(data)

        self.assertEqual(pickle_test.mean, 0.4)
        self.assertEqual(pickle_test.standard_dev, 0.1)

        self.assertEqual(pickle_test.lower_bound, 0)
        self.assertEqual(pickle_test.upper_bound, 1)

    def test_damping_sampling(self):
        test = msim.UncertainValue(2.2)
        test.set_distribution(msim.ParameterDistributionNormal(0, 1, 0.4, 0.1))
        test = msim.DampingSampling(test, 1, 2, 3, [0, 1], np.arange(2))

        data = pickle.dumps(test)
        pickle_test = pickle.loads(data)

        self.assertEqual(pickle_test.value.value, 2.2)

        self.assertEqual(pickle_test.value.get_distribution().mean, 0.4)
        self.assertEqual(
            pickle_test.value.get_distribution().standard_dev, 0.1)
        self.assertEqual(pickle_test.value.get_distribution().lower_bound, 0)
        self.assertEqual(pickle_test.value.get_distribution().upper_bound, 1)

        self.assertEqual(pickle_test.level, 1)
        self.assertEqual(pickle_test.type, 2)
        self.assertEqual(pickle_test.time, 3)

        self.assertEqual(pickle_test.matrix_indices, [0, 1])
        self.assertEqual(pickle_test.group_weights[0], 0)
        self.assertEqual(pickle_test.group_weights[1], 1)


if __name__ == '__main__':
    unittest.main()
