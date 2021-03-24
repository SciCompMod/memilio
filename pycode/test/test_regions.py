import unittest
import epidemiology.secir as secir

class Test_Regions(unittest.TestCase):
    def test_get_holidays(self):
        holidays = secir.get_holidays(9, start_date = secir.Date(2020, 10, 15), end_date = secir.Date(2020, 11, 15))
        self.assertEqual(len(holidays), 1)
        self.assertEqual(holidays[0][0], secir.Date(2020, 10, 31))
        self.assertEqual(holidays[0][1], secir.Date(2020, 11, 7))

    def test_state_id(self):
        self.assertEqual(secir.get_state_id(1001), 1)
        self.assertEqual(secir.get_state_id(1250), 9)

if __name__ == '__main__':
    unittest.main()