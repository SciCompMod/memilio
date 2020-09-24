import unittest
import epidemiology.secir as secir

class Test_UncertainValue(unittest.TestCase):
    def test_value(self):
        uv = secir.UncertainValue(0)
        self.assertEqual(uv.value, 0.0)
        uv.value = 1.0
        self.assertEqual(uv.value, 1.0)

    def test_distribution(self):
        uv = secir.UncertainValue(0)
        uv.set_distribution(secir.ParameterDistributionUniform(2, 2))
        uv.draw_sample()
        self.assertEqual(uv.value, 2.0)
        
if __name__ == '__main__':
    unittest.main()