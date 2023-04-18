"""
Module epidata_test provides tests for the epidata module.
"""

from memilio import progress_indicator

# Disable Progress Indicators for unittests
progress_indicator.ProgressIndicator.disable_indicators(True)
