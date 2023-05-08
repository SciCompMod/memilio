#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Patrick Lenz
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
import unittest
import time

from unittest.mock import patch
from io import StringIO

from memilio.epidata import progress_indicator

progress_indicator.ProgressIndicator.disable_indicators(True)


class Test_ProgressIndicator(unittest.TestCase):

    @patch('sys.stdout', new_callable=StringIO)
    def test_percentage_indicator(self, mock_print):
        progress_indicator.ProgressIndicator.disable_indicators(False)
        # test not full progress
        with progress_indicator.Percentage(delay=0) as p:
            p.set_progress(42/100)
        self.assertIn(' ]  42.00%', mock_print.getvalue())
        self.assertIn('[##', mock_print.getvalue())

        # test full progress
        with progress_indicator.Percentage(delay=0) as p:
            p.set_progress(100/100)
        self.assertIn('##] 100.00%', mock_print.getvalue())

        # test empty progress
        with progress_indicator.Percentage(delay=0) as p:
            p.set_progress(0/100)
        self.assertIn(' ]   0.00%', mock_print.getvalue())
        self.assertIn('[  ', mock_print.getvalue())

        progress_indicator.ProgressIndicator.disable_indicators(True)

    @patch('sys.stdout', new_callable=StringIO)
    def test_spinner(self, mock_print):
        s = progress_indicator.Spinner()
        spinner_prints = [" |", " /", " -", " \\"]
        for i in range(100):
            s.step()
            p = mock_print.getvalue()[-3:-1]
            self.assertEqual(p, spinner_prints[i % 4], 'help '+str(i))

    @patch('sys.stdout', new_callable=StringIO)
    def test_dots(self, mock_print):
        p = progress_indicator.Dots(message="testing dots", delay=0.1)
        for i in range(10):
            p.step()
            self.assertEqual(mock_print.getvalue(
            )[-17:-1], "testing dots ." + '.'*(i % 3) + " "*(2-(i % 3)))

    @patch('sys.stdout', new_callable=StringIO)
    def test_start_stop(self, mock_print):
        # progress indicator should be disabled
        self.assertEqual(progress_indicator.ProgressIndicator._disabled, True)
        # enable progress indicators
        progress_indicator.ProgressIndicator.disable_indicators(False)
        self.assertEqual(progress_indicator.ProgressIndicator._disabled, False)

        # test start and stop functionality
        p = progress_indicator.Dots(message="testing dots", delay=0.4)
        p.start()
        time.sleep(1.2)
        p.stop()
        output = mock_print.getvalue()
        # last output should have 3 dots
        self.assertEqual(output[-20:-4], "testing dots ...")

        # disable progress indicator for other tests
        progress_indicator.ProgressIndicator.disable_indicators(True)

    @patch('sys.stdout', new_callable=StringIO)
    def test_message_too_long(self, mock_print):
        p = progress_indicator.Percentage(delay=20, keep_output=True)
        p.start()
        p.set_message("A"*1000)
        p.stop()
        self.assertEqual(p._message, "")


if __name__ == '__main__':
    unittest.main()
