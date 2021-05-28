import unittest
import epidemiology.secir as secir
import numpy as np

class Test_UncertainMatrix(unittest.TestCase):
    def test_dampings(self):        
        m = secir.UncertainContactMatrix(secir.ContactMatrixGroup(num_matrices = 2, size = 2))
        d = secir.DampingSampling(
            value = secir.UncertainValue(3.0), 
            level = 1, 
            type = 2, 
            time = 10.0, 
            matrix_indices = [0, 1], 
            group_weights = np.r_[2.0, 1.0])
        m.dampings = [d]
        self.assertEqual(m.dampings[0].value.value, 3.0)
        self.assertEqual(m.dampings[0].level, 1)
        self.assertEqual(m.dampings[0].type, 2)
        self.assertEqual(m.dampings[0].time, 10.0)
        self.assertEqual(m.dampings[0].matrix_indices, [0, 1])
        self.assertTrue((m.dampings[0].group_weights == np.r_[2.0, 1.0]).all())
        
if __name__ == '__main__':
    unittest.main()