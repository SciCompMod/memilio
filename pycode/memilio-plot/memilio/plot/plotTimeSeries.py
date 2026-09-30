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
"""
:strong:`plotTimeSeries.py`
Standard plot of a MEmilio ``TimeSeries`` such as the result of a simulation.

The functions in this module accept a ``memilio.simulation.TimeSeries`` (or
any object with an ``as_ndarray()`` method returning an array of shape
``(1 + num_elements, num_time_points)`` whose first row holds the time
points) as well as a plain 2-D numpy array of that shape. The module does
not import ``memilio.simulation`` itself, so it can be used with the plot
package alone.
"""
from __future__ import annotations

import math
import warnings
from collections.abc import Sequence

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, NullFormatter

# Categorical palette in a fixed, colorblind-safe order. Compartment i is
# always drawn with _COLORS[i % 8]; from the ninth compartment on, the line
# style changes instead of introducing new hues.
_COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100',
           '#e87ba4', '#008300', '#4a3aa7', '#e34948']
_LINESTYLES = ['-', '--', '-.', ':']

_GRID_COLOR = '#e1e0d9'
_AXIS_COLOR = '#c3c2b7'
_TICK_COLOR = '#898781'
_TEXT_COLOR = '#52514e'
_TITLE_COLOR = '#0b0b0b'

# Legends with more rows than this are split into several columns.
_MAX_LEGEND_ROWS = 20


def _to_array(time_series) -> np.ndarray:
    """ Converts the input into a copied 2-D array (time points in row 0).

    :param time_series: TimeSeries-like object (with ``as_ndarray()``) or
        a 2-D array-like of shape (1 + num_elements, num_time_points).
    :returns: Copied numpy array of the same shape.
    """
    if hasattr(time_series, 'as_ndarray'):
        # as_ndarray() of the bindings is a view into C++ memory; copy it so
        # the plot data stays valid independent of the TimeSeries object.
        data = np.array(time_series.as_ndarray(), dtype=float, copy=True)
    else:
        data = np.array(time_series, dtype=float, copy=True)
    if data.ndim != 2 or data.shape[0] < 2:
        raise ValueError(
            'Expected a TimeSeries or a 2-D array of shape '
            '(1 + num_elements, num_time_points) with the time points in '
            f'the first row, got shape {data.shape}.')
    return data


def _names(items) -> list[str]:
    """ Converts a sequence of strings or objects with a ``name`` attribute
    (e.g. members of an InfectionState enum) into a list of strings.

    :param items: Iterable of strings or named objects, or a non-iterable
        object with a ``values()`` method returning such an iterable (for
        example the ``InfectionState`` enums of the bindings).
    :returns: List of names.
    """
    if isinstance(items, (str, bytes)):
        raise TypeError('Expected a sequence of names, not a single string.')
    if not hasattr(items, '__iter__'):
        if not callable(getattr(items, 'values', None)):
            raise TypeError(
                f'Expected a sequence of names, got {type(items).__name__}.')
        items = items.values()
    return [str(getattr(item, 'name', item)) for item in items]


def _prepare(time_series, labels, groups) -> tuple[
        np.ndarray, np.ndarray, list[str], list[str] | None]:
    """ Validates the input and reshapes the values by group and compartment.

    :param time_series: See :func:`plot_time_series`.
    :param labels: See :func:`plot_time_series`.
    :param groups: See :func:`plot_time_series`.
    :returns: Tuple of the time points (shape (num_time_points,)), the values
        (shape (num_groups, num_compartments, num_time_points)), the
        compartment names and the group names (None if no groups are given).
    """
    data = _to_array(time_series)
    times = data[0]
    values = data[1:]
    num_elements = values.shape[0]

    if groups is None:
        group_names = None
    elif isinstance(groups, (int, np.integer)) and not isinstance(
            groups, bool):
        if groups < 1:
            raise ValueError('groups must be a positive number of groups.')
        group_names = [f'Group {i}' for i in range(int(groups))]
    else:
        group_names = _names(groups)
        if not group_names:
            raise ValueError('groups must not be empty.')
    num_groups = 1 if group_names is None else len(group_names)

    if labels is None:
        if num_elements % num_groups != 0:
            raise ValueError(
                f'The number of elements ({num_elements}) is not divisible '
                f'by the number of groups ({num_groups}).')
        # Same default names as TimeSeries.print_table.
        labels = [f'C{i + 1}' for i in range(num_elements // num_groups)]
    else:
        labels = _names(labels)
        if len(labels) * num_groups != num_elements:
            if groups is None:
                raise ValueError(
                    f'Got {len(labels)} labels for {num_elements} elements. '
                    'Provide one label per element, or use the groups '
                    'argument if the elements are resolved by (age) groups.')
            raise ValueError(
                f'{len(labels)} labels times {num_groups} groups does not '
                f'match the number of elements ({num_elements}).')
    if len(set(labels)) != len(labels):
        raise ValueError('labels must be unique.')
    if group_names is not None and len(set(group_names)) != num_groups:
        raise ValueError('groups must be unique.')

    # Elements are ordered group-major: g0c0, g0c1, ..., g1c0, g1c1, ...
    values = values.reshape(num_groups, len(labels), values.shape[1])
    return times, values, labels, group_names


def _dates(times: np.ndarray, start_date) -> pd.DatetimeIndex:
    """ Converts simulation times (in days) into dates, rounded to seconds.

    :param times: Time points in days.
    :param start_date: Date (``datetime.date``, ``datetime.datetime``,
        ``pandas.Timestamp`` or ISO string) that corresponds to time 0.
    :returns: DatetimeIndex with one entry per time point.
    """
    dates = pd.Timestamp(start_date) + pd.to_timedelta(times, unit='D')
    return dates.round('s')


def time_series_to_dataframe(
        time_series, labels=None, groups=None, start_date=None
) -> pd.DataFrame:
    """ Converts a TimeSeries into a tidy (long-form) pandas DataFrame.

    The frame has one row per time point and element with the columns
    ``Time``, ``Date`` (only if ``start_date`` is given), ``Groups`` (only
    if ``groups`` is given), ``Compartments`` and ``Values``.
    ``Compartments`` and ``Groups`` are ordered categoricals in the order of
    the TimeSeries, so the frame can be used directly with
    grammar-of-graphics libraries such as seaborn, plotnine or altair. A wide
    table (one column per compartment, summed over groups) is obtained by
    ``df.pivot_table(index='Time', columns='Compartments', values='Values',
    aggfunc='sum')``.

    :param time_series: ``memilio.simulation.TimeSeries`` (or any object with
        ``as_ndarray()``) or a 2-D array of shape
        (1 + num_elements, num_time_points) with the time points in row 0.
    :param labels: Names of the compartments (elements). Either strings or
        objects with a ``name`` attribute, e.g. ``InfectionState.values()``
        of a model. If None, the elements are named ``C1``, ``C2``, ...
        (Default value = None)
    :param groups: Number of groups or their names (e.g. age groups) if the
        elements are resolved by groups. The elements are expected in the
        order of the bindings, i.e. all compartments of the first group,
        then all compartments of the second group, and so on.
        (Default value = None)
    :param start_date: Date corresponding to time 0. If given, a ``Date``
        column (rounded to full seconds) is added. (Default value = None)
    :returns: Long-form DataFrame.
    """
    times, values, labels, group_names = _prepare(
        time_series, labels, groups)
    num_groups, num_compartments, num_time_points = values.shape

    frame = {'Time': np.tile(times, num_groups * num_compartments)}
    if start_date is not None:
        frame['Date'] = _dates(frame['Time'], start_date)
    if group_names is not None:
        frame['Groups'] = pd.Categorical.from_codes(
            np.repeat(np.arange(num_groups),
                      num_compartments * num_time_points),
            categories=group_names, ordered=True)
    frame['Compartments'] = pd.Categorical.from_codes(
        np.tile(np.repeat(np.arange(num_compartments), num_time_points),
                num_groups),
        categories=labels, ordered=True)
    frame['Values'] = values.reshape(-1)
    return pd.DataFrame(frame)


def _resolve_select(select, labels: list[str]) -> list[int]:
    """ Resolves the compartments to plot into sorted, unique indices.

    :param select: None (all), a single name/index or a sequence of names
        or indices.
    :param labels: All compartment names.
    :returns: Sorted list of compartment indices.
    """
    if select is None:
        return list(range(len(labels)))
    if isinstance(select, (str, int, np.integer)) or hasattr(select, 'name'):
        select = [select]
    indices = set()
    for item in select:
        if isinstance(item, (int, np.integer)) and not isinstance(item, bool):
            if not 0 <= item < len(labels):
                raise ValueError(f'Compartment index {item} out of range.')
            indices.add(int(item))
        else:
            name = str(getattr(item, 'name', item))
            if name not in labels:
                raise ValueError(
                    f'Unknown compartment {name!r}. '
                    f'Available: {", ".join(labels)}.')
            indices.add(labels.index(name))
    return sorted(indices)


def _color_for(index: int, name: str, colors) -> str:
    """ Color of a line, keyed by the position of its compartment (or group).

    :param index: Index of the compartment (or group).
    :param name: Name of the compartment (or group).
    :param colors: None, a sequence of colors (by index) or a dictionary
        mapping names to colors.
    :returns: Matplotlib color.
    """
    if isinstance(colors, dict):
        if name in colors:
            return colors[name]
    elif colors is not None and index < len(colors):
        return colors[index]
    return _COLORS[index % len(_COLORS)]


def _style_for(index: int) -> str:
    """ Line style of a line, changing once the palette has been used up.

    :param index: Index of the compartment (or group).
    :returns: Matplotlib line style.
    """
    return _LINESTYLES[(index // len(_COLORS)) % len(_LINESTYLES)]


def _format_count(value, _position=None) -> str:
    """ Tick formatter with thousands separators for large numbers.

    :param value: Tick value.
    :param _position: Tick position (unused). (Default value = None)
    :returns: Formatted tick label.
    """
    if abs(value) >= 1000:
        return f'{value:,.0f}'
    if 0 < abs(value) < 1e-3:
        return f'{value:.0e}'.replace('e-0', 'e-')
    return f'{value:g}'


def _apply_style(ax, log_scale: bool, use_dates: bool):
    """ Applies the style to an axes: recessive grid and spines,
    readable tick labels.

    :param ax: Matplotlib axes.
    :param log_scale: Whether the y axis is logarithmic.
    :param use_dates: Whether the x axis shows dates.
    """
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(_AXIS_COLOR)
    ax.tick_params(colors=_TICK_COLOR, labelcolor=_TEXT_COLOR)
    ax.xaxis.label.set_color(_TEXT_COLOR)
    ax.yaxis.label.set_color(_TEXT_COLOR)
    ax.grid(True, axis='y', color=_GRID_COLOR, linewidth=0.8, linestyle='-')
    ax.grid(False, axis='x')
    ax.set_axisbelow(True)
    ax.margins(x=0)

    if log_scale:
        ax.set_yscale('log', nonpositive='clip')
        ax.yaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(FuncFormatter(_format_count))

    if use_dates:
        locator = mdates.AutoDateLocator()
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))


def _draw_legend(ax, num_lines: int):
    """ Draws the legend right of the axes, split into columns if it has
    many entries.

    :param ax: Matplotlib axes.
    :param num_lines: Number of legend entries.
    """
    ncol = max(1, math.ceil(num_lines / _MAX_LEGEND_ROWS))
    ax.legend(frameon=False, ncol=ncol, loc='upper left',
              bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)


def _collect_lines(values: np.ndarray, labels: list[str],
                   group_names: list[str] | None, selected: list[int],
                   sum_groups: bool, colors) -> list[tuple]:
    """ Determines the lines to draw.

    :param values: Values of shape (num_groups, num_compartments,
        num_time_points).
    :param labels: Compartment names.
    :param group_names: Group names or None.
    :param selected: Indices of the compartments to draw.
    :param sum_groups: Whether to sum the compartments over the groups.
    :param colors: See :func:`plot_time_series`.
    :returns: List of (color, line style, label, y values) per line.
    """
    lines = []
    for index in selected:
        name = labels[index]
        if group_names is None or sum_groups:
            lines.append((_color_for(index, name, colors), _style_for(index),
                          name, values[:, index, :].sum(axis=0)))
        elif len(selected) == 1:
            # A single compartment: the groups take the colors.
            for group_index, group_name in enumerate(group_names):
                lines.append((_color_for(group_index, group_name, colors),
                              _style_for(group_index),
                              f'{name} ({group_name})',
                              values[group_index, index, :]))
        else:
            for group_index, group_name in enumerate(group_names):
                lines.append((_color_for(index, name, colors),
                              _LINESTYLES[group_index % len(_LINESTYLES)],
                              f'{name} ({group_name})',
                              values[group_index, index, :]))
    styles = [(color, style) for color, style, _, _ in lines]
    if len(set(styles)) < len(styles):
        warnings.warn(
            'Some lines share color and line style and cannot be told '
            'apart. Use select to plot fewer compartments or groups.',
            stacklevel=3)
    return lines


def plot_time_series(
        time_series, labels=None, *, groups=None, sum_groups: bool = True,
        select=None, ax=None, title: str | None = None,
        xlabel: str | None = None, ylabel: str = 'Number of individuals',
        log_scale: bool = False, start_date=None,
        legend: bool = True,
        colors: Sequence[str] | dict[str, str] | None = None,
        figsize: tuple[float, float] = (8, 4.5), **plot_kwargs):
    """ Plots the elements of a TimeSeries over time, one line per compartment.

    Example::

        from memilio.simulation.oseir import InfectionState, simulate
        from memilio.plot.plotTimeSeries import plot_time_series

        # model: oseir.Model with two age groups
        result = simulate(0, 100, 0.1, model)
        ax = plot_time_series(result, labels=InfectionState.values(),
                              groups=['0-19', '20+'], title='ODE SEIR')
        ax.figure.savefig('seir.pdf')

    Every compartment has a fixed color. For more than eight compartments,
    the colors are reused with a different line style; consider ``select``
    to plot only the compartments of interest.

    :param time_series: ``memilio.simulation.TimeSeries`` (or any object with
        ``as_ndarray()``) or a 2-D array of shape
        (1 + num_elements, num_time_points) with the time points in row 0.
    :param labels: Names of the compartments (elements). Either strings or
        objects with a ``name`` attribute, e.g. ``InfectionState.values()``
        of a model. If None, the elements are named ``C1``, ``C2``, ...
        (Default value = None)
    :param groups: Number of groups or their names (e.g. age groups) if the
        elements are resolved by groups, i.e. if the TimeSeries has
        ``len(labels) * num_groups`` elements ordered group by group.
        (Default value = None)
    :param sum_groups: If True, the compartments are summed over all groups
        and one line per compartment is drawn. If False, one line per group
        and compartment is drawn, labeled ``'<compartment> (<group>)'``. The
        groups are distinguished by line style, or by color if a single
        compartment is selected. A warning is issued if the lines cannot be
        told apart. (Default value = True)
    :param select: Compartment(s) to plot, given by name or index. None plots
        all compartments. Colors are assigned by the position of a compartment
        in the TimeSeries, so a selection does not change the colors.
        (Default value = None)
    :param ax: Matplotlib axes to draw into. If None, a new figure with
        constrained layout is created. When drawing into your own figure,
        create it with ``layout='constrained'`` or save it with
        ``bbox_inches='tight'`` so that the legend right of the axes is not
        cut off. (Default value = None)
    :param title: Title of the plot. (Default value = None)
    :param xlabel: Label of the x axis. If None, 'Time (days)' or 'Date' is
        used. (Default value = None)
    :param ylabel: Label of the y axis.
        (Default value = 'Number of individuals')
    :param log_scale: Use a logarithmic y axis. (Default value = False)
    :param start_date: Date corresponding to time 0 (``datetime.date`` or
        similar). If given, the x axis shows dates. (Default value = None)
    :param legend: Whether to draw a legend. It is placed right of the axes.
        (Default value = True)
    :param colors: Sequence of colors indexed by compartment, or a dictionary
        mapping compartment names to colors. Compartments not covered use the
        default palette. If a single compartment is plotted per group, the
        colors refer to the groups instead. (Default value = None)
    :param figsize: Size of the figure in inches if a new figure is created.
        (Default value = (8, 4.5))
    :param plot_kwargs: Additional keyword arguments passed to
        ``matplotlib.axes.Axes.plot`` for every line, e.g. ``linewidth``.
    :returns: The matplotlib axes containing the plot. Use ``ax.figure`` to
        access and save the figure.
    """
    times, values, labels, group_names = _prepare(
        time_series, labels, groups)
    selected = _resolve_select(select, labels)
    x = _dates(times, start_date) if start_date is not None else times

    lines = _collect_lines(values, labels, group_names, selected,
                           sum_groups, colors)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize, layout='constrained')

    line_kwargs = {'linewidth': 2.0, 'solid_capstyle': 'round',
                   'solid_joinstyle': 'round'}
    for alias, full_name in (('lw', 'linewidth'), ('ls', 'linestyle'),
                             ('c', 'color')):
        if alias in plot_kwargs:
            plot_kwargs[full_name] = plot_kwargs.pop(alias)
    plot_kwargs.pop('label', None)
    line_kwargs.update(plot_kwargs)
    for color, style, name, y in lines:
        ax.plot(x, y, label=name,
                **{'color': color, 'linestyle': style, **line_kwargs})

    _apply_style(ax, log_scale, start_date is not None)

    if xlabel is None:
        xlabel = 'Date' if start_date is not None else 'Time (days)'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title is not None:
        ax.set_title(title, fontweight='bold', color=_TITLE_COLOR)
    if legend and lines:
        _draw_legend(ax, len(lines))
    return ax
