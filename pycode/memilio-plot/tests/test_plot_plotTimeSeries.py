#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Kilian Volmer
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
import datetime as dt
import enum
import unittest
import warnings

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.layout_engine import ConstrainedLayoutEngine
from matplotlib.ticker import NullFormatter

import memilio.plot.plotTimeSeries as pts

try:
    import memilio.simulation as mio
except ImportError:
    # The plot package is tested without the compiled bindings.
    mio = None

matplotlib.use('Agg')


class _StubTimeSeries:
    """Minimal stand-in for memilio.simulation.TimeSeries."""

    def __init__(self, data):
        self._data = np.asarray(data, dtype=float)

    def as_ndarray(self):
        """Returns the data as array (1 + num_elements, num_time_points)."""
        return self._data


class _State(enum.Enum):
    """Iterable enum with named members, like Python's own enums."""
    Susceptible = 0
    Infected = 1
    Recovered = 2


class _NamedMember:
    """Object with a name attribute, like a member of a bindings enum."""

    def __init__(self, name):
        self.name = name


class _BindingsLikeEnum:
    """Mimics the InfectionState enums of the bindings: not iterable, but
    with a static values() method returning members with a name."""

    @staticmethod
    def values():
        """Returns the members of the enum."""
        return [_NamedMember('Susceptible'), _NamedMember('Infected'),
                _NamedMember('Recovered')]


class TestPlotTimeSeries(unittest.TestCase):
    """Tests for plot_time_series and time_series_to_dataframe."""

    def setUp(self):
        """Creates a TimeSeries with 2 groups x 3 compartments = 6 elements
        in group-major order, where element i has the values (i+1)*10 + t."""
        self.times = np.array([0.0, 0.5, 1.0, 2.0, 3.0])
        self.num_compartments = 3
        self.values = np.array(
            [(i + 1) * 10.0 + self.times for i in range(6)])
        self.data = np.vstack([self.times, self.values])
        self.ts = _StubTimeSeries(self.data)
        self.labels = ['S', 'I', 'R']
        self.groups = ['0-4', '5+']

    def tearDown(self):
        """Closes all figures created by the test."""
        plt.close('all')

    def _labels(self, ax):
        """Returns the labels of the lines in drawing order."""
        return [line.get_label() for line in ax.get_lines()]

    def _legend_texts(self, ax):
        """Returns the entries of the legend."""
        return [t.get_text() for t in ax.get_legend().get_texts()]

    def _legend_is_outside(self, ax):
        """Returns whether the legend is anchored right of the axes."""
        ax.figure.canvas.draw()
        anchor = ax.get_legend().get_bbox_to_anchor()
        return anchor.x0 >= ax.get_window_extent().x1

    # ---- input handling -------------------------------------------------

    def test_default_labels(self):
        """Test that elements without labels are named C1, C2, ... as in
        TimeSeries.print_table, and that the legend lists them."""
        ax = pts.plot_time_series(self.ts)
        expected = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6']
        self.assertEqual(self._labels(ax), expected)
        self.assertEqual(self._legend_texts(ax), expected)

    def test_given_labels(self):
        """Test that one line per element is drawn with the given labels and
        the unchanged time points and values, for a TimeSeries and for a
        plain array."""
        labels = ['a', 'b', 'c', 'd', 'e', 'f']
        for time_series in (self.ts, self.data):
            ax = pts.plot_time_series(time_series, labels)
            self.assertEqual(self._labels(ax), labels)
            for line, row in zip(ax.get_lines(), self.values):
                np.testing.assert_array_equal(line.get_xdata(), self.times)
                np.testing.assert_array_equal(line.get_ydata(), row)

    def test_labels_from_enum(self):
        """Test that labels can be taken from the members of an enum, from
        a bindings-like enum class or from its values()."""
        data = self.data[:4]  # 3 elements
        expected = ['Susceptible', 'Infected', 'Recovered']
        for labels in (_State, _BindingsLikeEnum, _BindingsLikeEnum.values()):
            ax = pts.plot_time_series(_StubTimeSeries(data), labels)
            self.assertEqual(self._labels(ax), expected)

    def test_input_is_copied(self):
        """Test that the plot keeps its data if the input array is modified
        afterwards, since as_ndarray() of the bindings is a view."""
        data = self.data.copy()
        ax = pts.plot_time_series(_StubTimeSeries(data), groups=2)
        data[1:] = 0.0
        self.assertGreater(ax.get_lines()[0].get_ydata().sum(), 0.0)

    def test_invalid_input(self):
        """Test that wrong shapes, mismatching numbers of labels and groups,
        duplicate names and wrong argument types raise errors."""
        cases = [
            ('1-D input', {'time_series': np.array([0.0, 1.0, 2.0])},
             ValueError),
            ('no elements', {'time_series': np.zeros((1, 3))}, ValueError),
            ('too few labels', {'labels': ['only', 'two']}, ValueError),
            ('labels times groups mismatch',
             {'labels': self.labels, 'groups': 3}, ValueError),
            ('group names mismatch',
             {'labels': self.labels, 'groups': ['x']}, ValueError),
            ('elements not divisible by groups', {'groups': 4}, ValueError),
            ('zero groups', {'groups': 0}, ValueError),
            ('empty groups', {'groups': []}, ValueError),
            ('duplicate labels',
             {'labels': ['a', 'a', 'b', 'c', 'd', 'e']}, ValueError),
            ('labels as single string', {'labels': 'abcdef'}, TypeError),
            ('groups as bool', {'groups': True}, TypeError),
            ('groups as float', {'groups': 2.0}, TypeError),
        ]
        for description, kwargs, error in cases:
            with self.subTest(description):
                kwargs = {'time_series': self.ts, **kwargs}
                with self.assertRaises(error):
                    pts.plot_time_series(**kwargs)

    # ---- groups ---------------------------------------------------------

    def test_sum_groups(self):
        """Test that the compartments are summed over the groups by
        default."""
        ax = pts.plot_time_series(self.ts, self.labels, groups=self.groups)
        self.assertEqual(self._labels(ax), self.labels)
        for c, line in enumerate(ax.get_lines()):
            expected = self.values[c] + self.values[c + self.num_compartments]
            np.testing.assert_allclose(line.get_ydata(), expected)

    def test_groups_by_count(self):
        """Test that groups given as a number are named Group 0, Group 1,
        ... and that default labels are per compartment."""
        ax = pts.plot_time_series(self.ts, groups=2)
        self.assertEqual(self._labels(ax), ['C1', 'C2', 'C3'])
        ax = pts.plot_time_series(self.ts, groups=2, sum_groups=False)
        self.assertEqual(self._labels(ax)[:2],
                         ['C1 (Group 0)', 'C1 (Group 1)'])

    def test_per_group_lines(self):
        """Test one line per group and compartment: the color follows the
        compartment and the line style the group; without groups,
        sum_groups has no effect."""
        ax = pts.plot_time_series(
            self.ts, self.labels, groups=self.groups, sum_groups=False)
        lines = ax.get_lines()
        self.assertEqual(
            self._labels(ax),
            ['S (0-4)', 'S (5+)', 'I (0-4)', 'I (5+)', 'R (0-4)', 'R (5+)'])
        self.assertEqual(lines[0].get_color(), lines[1].get_color())
        self.assertNotEqual(lines[0].get_color(), lines[2].get_color())
        self.assertNotEqual(lines[0].get_linestyle(), lines[1].get_linestyle())
        self.assertEqual(lines[0].get_linestyle(), lines[2].get_linestyle())
        np.testing.assert_allclose(lines[1].get_ydata(), self.values[3])

        ax = pts.plot_time_series(self.ts, sum_groups=False)
        self.assertEqual(len(ax.get_lines()), 6)

    def test_single_compartment_per_group(self):
        """Test that the groups take the colors (default or given) if a
        single compartment is plotted per group."""
        ax = pts.plot_time_series(
            self.ts, self.labels, groups=self.groups, sum_groups=False,
            select='I')
        lines = ax.get_lines()
        self.assertEqual(self._labels(ax), ['I (0-4)', 'I (5+)'])
        self.assertEqual(lines[0].get_color(), '#2a78d6')
        self.assertEqual(lines[1].get_color(), '#eb6834')
        self.assertEqual(lines[0].get_linestyle(), lines[1].get_linestyle())
        np.testing.assert_allclose(lines[0].get_ydata(), self.values[1])
        np.testing.assert_allclose(lines[1].get_ydata(), self.values[4])

        ax = pts.plot_time_series(
            self.ts, self.labels, groups=self.groups, sum_groups=False,
            select='I', colors=['black', 'gray'])
        self.assertEqual([line.get_color() for line in ax.get_lines()],
                         ['black', 'gray'])

    def test_ambiguous_lines_warn(self):
        """Test that a warning is issued if lines share color and line
        style (more groups than line styles), and not otherwise."""
        data = np.vstack([self.times, np.ones((10, len(self.times)))])
        with self.assertWarns(UserWarning):
            pts.plot_time_series(data, ['a', 'b'], groups=5, sum_groups=False)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            pts.plot_time_series(data, ['a', 'b'], groups=5)
            pts.plot_time_series(
                data, ['a', 'b'], groups=5, sum_groups=False, select='a')

    # ---- selection and colors ------------------------------------------

    def test_select(self):
        """Test selecting compartments by name, index or enum member, in the
        order of the TimeSeries and without changing their colors."""
        ax_all = pts.plot_time_series(self.ts, self.labels, groups=2)
        colors = {line.get_label(): line.get_color()
                  for line in ax_all.get_lines()}
        ax = pts.plot_time_series(
            self.ts, self.labels, groups=2, select=['R', 'I'])
        self.assertEqual(self._labels(ax), ['I', 'R'])
        for line in ax.get_lines():
            self.assertEqual(line.get_color(), colors[line.get_label()])
        ax = pts.plot_time_series(self.ts, self.labels, groups=2, select=2)
        self.assertEqual(self._labels(ax), ['R'])
        ax = pts.plot_time_series(
            self.ts, self.labels, groups=2, select=_NamedMember('I'))
        self.assertEqual(self._labels(ax), ['I'])
        ax = pts.plot_time_series(
            self.ts, self.labels, groups=2, select=[0, 'S', np.int64(1)])
        self.assertEqual(self._labels(ax), ['S', 'I'])
        with self.assertRaises(ValueError):
            pts.plot_time_series(self.ts, self.labels, groups=2, select='X')
        with self.assertRaises(ValueError):
            pts.plot_time_series(self.ts, self.labels, groups=2, select=[7])

    def test_custom_colors(self):
        """Test overriding the colors by a dictionary of compartment names
        or by a sequence indexed by compartment."""
        ax = pts.plot_time_series(
            self.ts, self.labels, groups=2, colors={'I': 'black'})
        self.assertEqual([line.get_color() for line in ax.get_lines()],
                         ['#2a78d6', 'black', '#1baf7a'])
        ax = pts.plot_time_series(
            self.ts, self.labels, groups=2, colors=['red', 'green'])
        self.assertEqual([line.get_color() for line in ax.get_lines()],
                         ['red', 'green', '#1baf7a'])

    def test_more_than_eight_compartments(self):
        """Test that from the ninth compartment on, the colors are reused
        with a different line style."""
        data = np.vstack([self.times, np.ones((10, len(self.times)))])
        ax = pts.plot_time_series(data)
        lines = ax.get_lines()
        self.assertEqual(len(lines), 10)
        self.assertEqual(lines[8].get_color(), lines[0].get_color())
        self.assertEqual(lines[0].get_linestyle(), '-')
        self.assertEqual(lines[8].get_linestyle(), '--')
        self.assertEqual(len(self._legend_texts(ax)), 10)

    # ---- legend ---------------------------------------------------------

    def test_legend_right_of_axes(self):
        """Test that the legend is placed right of the axes, also for a
        single line and for a user-provided axes."""
        ax = pts.plot_time_series(self.ts, groups=2)
        self.assertTrue(self._legend_is_outside(ax))
        ax = pts.plot_time_series(self.ts, self.labels, groups=2, select='I')
        self.assertEqual(self._legend_texts(ax), ['I'])
        self.assertTrue(self._legend_is_outside(ax))
        _, axes = plt.subplots(2, 2)
        ax = pts.plot_time_series(self.ts, ax=axes[0, 0])
        self.assertEqual(len(self._legend_texts(ax)), 6)
        self.assertTrue(self._legend_is_outside(ax))

    def test_no_legend(self):
        """Test that no legend is drawn if it is disabled or if there are no
        lines."""
        ax = pts.plot_time_series(self.ts, groups=2, legend=False)
        self.assertIsNone(ax.get_legend())
        ax = pts.plot_time_series(self.ts, groups=2, select=[])
        self.assertEqual(len(ax.get_lines()), 0)
        self.assertIsNone(ax.get_legend())

    # ---- axes, labels and options --------------------------------------

    def test_new_figure(self):
        """Test that a new figure of the given size with constrained layout
        and default axis labels is created if no axes is given."""
        ax = pts.plot_time_series(self.ts, groups=2, figsize=(10, 3))
        np.testing.assert_array_equal(
            ax.figure.get_size_inches(), [10.0, 3.0])
        self.assertIsInstance(
            ax.figure.get_layout_engine(), ConstrainedLayoutEngine)
        self.assertEqual(ax.get_xlabel(), 'Time (days)')
        self.assertEqual(ax.get_ylabel(), 'Number of individuals')
        self.assertEqual(ax.get_title(), '')

    def test_given_axes_and_labels(self):
        """Test drawing into a given axes with a title, axis labels and line
        keyword arguments."""
        _, ax_in = plt.subplots()
        ax = pts.plot_time_series(
            self.ts, self.labels, groups=2, ax=ax_in, title='T',
            xlabel='x', ylabel='y', linewidth=4.0)
        self.assertIs(ax, ax_in)
        self.assertEqual(ax.get_title(), 'T')
        self.assertEqual(ax.get_xlabel(), 'x')
        self.assertEqual(ax.get_ylabel(), 'y')
        self.assertEqual(ax.get_yscale(), 'linear')
        self.assertFalse(ax.spines['top'].get_visible())
        for line in ax.get_lines():
            self.assertEqual(line.get_linewidth(), 4.0)

    def test_plot_kwargs(self):
        """Test that keyword arguments and their matplotlib aliases are
        passed to all lines and override the defaults, except the label."""
        ax = pts.plot_time_series(
            self.ts, groups=2, lw=3.0, color='black', ls=':', label='x')
        for line in ax.get_lines():
            self.assertEqual(line.get_linewidth(), 3.0)
            self.assertEqual(line.get_color(), 'black')
            self.assertEqual(line.get_linestyle(), ':')
        self.assertEqual(self._labels(ax), ['C1', 'C2', 'C3'])

    def test_log_scale(self):
        """Test the logarithmic y axis with the count formatter on the major
        ticks and no labels on the minor ticks."""
        ax = pts.plot_time_series(self.ts, groups=2, log_scale=True)
        self.assertEqual(ax.get_yscale(), 'log')
        self.assertEqual(ax.yaxis.get_major_formatter()(1000.0), '1,000')
        self.assertIsInstance(ax.yaxis.get_minor_formatter(), NullFormatter)

    def test_start_date(self):
        """Test that the time points are converted to dates from a date
        object or an ISO string, with 'Date' as axis label."""
        for start_date in (dt.date(2020, 3, 1), '2020-03-01'):
            ax = pts.plot_time_series(
                self.ts, groups=2, start_date=start_date)
            self.assertEqual(ax.get_xlabel(), 'Date')
            xdata = pd.DatetimeIndex(ax.get_lines()[0].get_xdata())
            self.assertEqual(xdata[0], pd.Timestamp('2020-03-01'))
            self.assertEqual(xdata[1], pd.Timestamp('2020-03-01 12:00'))
            self.assertEqual(xdata[-1], pd.Timestamp('2020-03-04'))

    def test_tick_formatter(self):
        """Test the y tick labels: thousands separators for large numbers,
        plain fractions and a consistent notation for very small values."""
        ax = pts.plot_time_series(self.ts, groups=2)
        formatter = ax.yaxis.get_major_formatter()
        self.assertEqual(formatter(1234567.0), '1,234,567')
        self.assertEqual(formatter(12.5), '12.5')
        self.assertEqual(formatter(0.0), '0')
        self.assertEqual(formatter(0.001), '0.001')
        self.assertEqual(formatter(1e-4), '1e-4')
        self.assertEqual(formatter(1e-5), '1e-5')

    # ---- dataframe ------------------------------------------------------

    def test_to_dataframe(self):
        """Test the long-form data frame with groups: columns, categorical
        order, values per group and compartment, and summing by pivoting."""
        df = pts.time_series_to_dataframe(self.ts, self.labels, self.groups)
        self.assertEqual(list(df.columns),
                         ['Time', 'Groups', 'Compartments', 'Values'])
        self.assertEqual(len(df), 6 * len(self.times))
        self.assertEqual(list(df['Compartments'].cat.categories), self.labels)
        self.assertEqual(list(df['Groups'].cat.categories), self.groups)
        self.assertTrue(df['Compartments'].cat.ordered)
        # Element 4 = second group, second compartment.
        subset = df[(df['Groups'] == '5+') & (df['Compartments'] == 'I')]
        np.testing.assert_array_equal(subset['Time'].to_numpy(), self.times)
        np.testing.assert_allclose(
            subset['Values'].to_numpy(), self.values[4])
        wide = df.pivot_table(index='Time', columns='Compartments',
                              values='Values', aggfunc='sum', observed=False)
        np.testing.assert_allclose(
            wide['S'].to_numpy(), self.values[0] + self.values[3])

    def test_to_dataframe_without_groups_with_dates(self):
        """Test the data frame without groups and with a date column that is
        rounded to full seconds."""
        data = np.vstack([[0.0, 1.253706], np.ones((2, 2))])
        df = pts.time_series_to_dataframe(
            data, start_date=dt.date(2021, 1, 1))
        self.assertEqual(list(df.columns),
                         ['Time', 'Date', 'Compartments', 'Values'])
        self.assertEqual(list(df['Compartments'].cat.categories),
                         ['C1', 'C2'])
        self.assertEqual(df['Date'].iloc[0], pd.Timestamp('2021-01-01'))
        self.assertEqual(df['Date'].iloc[1],
                         pd.Timestamp('2021-01-02 06:05:20'))

    # ---- real bindings --------------------------------------------------

    @unittest.skipIf(mio is None, 'memilio.simulation not installed')
    def test_with_bindings_time_series(self):
        """Test plotting and converting a TimeSeries of the bindings; skipped
        if memilio.simulation is not installed."""
        ts = mio.TimeSeries(4)
        ts.add_time_point(0.0, np.r_[100.0, 10.0, 200.0, 20.0])
        ts.add_time_point(1.0, np.r_[90.0, 20.0, 180.0, 40.0])
        ax = pts.plot_time_series(ts, ['S', 'I'], groups=['a', 'b'])
        lines = ax.get_lines()
        self.assertEqual(self._labels(ax), ['S', 'I'])
        np.testing.assert_allclose(lines[0].get_ydata(), [300.0, 270.0])
        np.testing.assert_allclose(lines[1].get_ydata(), [30.0, 60.0])
        df = pts.time_series_to_dataframe(ts, ['S', 'I'], groups=2)
        self.assertEqual(len(df), 8)


if __name__ == '__main__':
    unittest.main()
