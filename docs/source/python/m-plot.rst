MEmilio Plot
=============

MEmilio Plot provides modules and scripts to plot epidemiological or simulation data as returned
by other packages of the MEmilio software.

The package is contained inside the folder `pycode/memilio-plot <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-plot>`_.

.. note:: This package is under active development.

Installation
------------

See :ref:`python-package-installation` for a detailed installation guide.

Dependencies
------------

Required python packages:

- pandas>=1.2.2
- matplotlib>=3.6
- numpy>=1.22, !=1.25.*
- openpyxl
- xlrd
- requests
- pyxlsb
- wget
- folium
- mapclassify
- geopandas
- h5py
- imageio
- datetime

Plotting simulation results
---------------------------

The module ``memilio.plot.plotTimeSeries`` provides a standard plot for ``TimeSeries`` objects as returned by the
simulations of the :doc:`MEmilio Python bindings <m-simulation>`:

.. code-block:: python

    from datetime import date

    import numpy as np

    from memilio.plot.plotTimeSeries import plot_time_series
    from memilio.simulation import AgeGroup, Damping
    from memilio.simulation.osecir import InfectionState, Model, simulate

    # ODE SECIR model with two age groups
    groups = ['0-19', '20+']
    populations = [15000, 68000]
    num_groups = len(groups)

    model = Model(num_groups)
    for i, population in enumerate(populations):
        group = AgeGroup(i)
        # time spent in the compartments (days)
        model.parameters.TimeExposed[group] = 3.2
        model.parameters.TimeInfectedNoSymptoms[group] = 2.
        model.parameters.TimeInfectedSymptoms[group] = 6.
        model.parameters.TimeInfectedSevere[group] = 12.
        model.parameters.TimeInfectedCritical[group] = 8.
        # transmission and transition probabilities
        model.parameters.TransmissionProbabilityOnContact[group] = 1.0
        model.parameters.RelativeTransmissionNoSymptoms[group] = 0.67
        model.parameters.RiskOfInfectionFromSymptomatic[group] = 0.25
        model.parameters.MaxRiskOfInfectionFromSymptomatic[group] = 0.5
        model.parameters.RecoveredPerInfectedNoSymptoms[group] = 0.09
        model.parameters.SeverePerInfectedSymptoms[group] = 0.2
        model.parameters.CriticalPerSevere[group] = 0.25
        model.parameters.DeathsPerCritical[group] = 0.3
        # initial populations
        model.populations[group, InfectionState.Exposed] = 100
        model.populations[group, InfectionState.InfectedNoSymptoms] = 50
        model.populations[group, InfectionState.InfectedSymptoms] = 50
        model.populations[group, InfectionState.InfectedSevere] = 20
        model.populations[group, InfectionState.InfectedCritical] = 10
        model.populations[group, InfectionState.Recovered] = 10
        model.populations.set_difference_from_group_total_AgeGroup(
            (group, InfectionState.Susceptible), population)

    # contacts: one contact per person and day, reduced by 90% from day 30 on
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (num_groups, num_groups))
    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=np.ones((num_groups, num_groups)) * 0.9, t=30., level=0, type=0))
    model.check_constraints()

    result = simulate(0, 100, 0.1, model)
    ax = plot_time_series(
        result, labels=InfectionState.values(), groups=groups,
        select=['Exposed', 'InfectedSymptoms', 'Dead'],
        start_date=date(2020, 3, 1), title='ODE SECIR simulation')
    ax.figure.savefig('secir.pdf')

- ``labels`` names the compartments. It accepts strings or the ``InfectionState.values()`` of a model. Without labels,
  the elements are named ``C1``, ``C2``, ... .
- ``groups`` (the number or the names of, e.g., age groups) tells the function how the elements of the ``TimeSeries``
  are arranged. By default, the compartments are summed over all groups; with ``sum_groups=False`` one line per group
  and compartment is drawn.
- ``select`` restricts the plot to some compartments, ``log_scale`` switches to a logarithmic axis and ``start_date``
  turns the time axis into calendar dates. Colors are fixed per compartment, so a selection does not change them.
- The function returns a matplotlib ``Axes``. Pass ``ax=`` to draw into an existing figure (e.g., a subplot grid) and
  use ``ax.figure.savefig`` to save the figure. The legend is placed right of the axes, so create your own figures with
  ``layout='constrained'`` or save them with ``bbox_inches='tight'``.

``time_series_to_dataframe`` converts a ``TimeSeries`` into a tidy pandas ``DataFrame`` with the columns ``Time``,
``Groups``, ``Compartments`` and ``Values`` (and ``Date`` if a start date is given) for use with other libraries such as
seaborn, plotnine or altair. A complete example is given in
`examples/plot/plotSimulationResults.py <https://github.com/SciCompMod/memilio/blob/main/pycode/examples/plot/plotSimulationResults.py>`_.
