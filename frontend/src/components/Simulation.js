import React, { PureComponent } from 'react';
import { connect } from 'react-redux';
import { Button } from 'reactstrap';
import { getActiveMeasures } from '../redux/measures';
import { getParameterMap } from '../redux/parameters';
import {getPopulationsOfRegion, getSelectedChildData, getSelectedData} from '../redux/app';
import {setData, setRegionData} from '../redux/seir';

import { simulate_seir, makeSeirParam, Damping } from '../common/seir.js';
import { calculateDamping } from '../common/utils';

class Simulation extends PureComponent {
  simulate() {
    const selectedResult = this.simulateRegion(this.props.start, this.props.selected.population);
    this.props.setData(selectedResult);

    const regionResults = {};
    for (let [regionId, region] of Object.entries(this.props.childData.all)) {
      const start = region[0]
      const population = this.props.populations[regionId].EWZ;
      regionResults[regionId] = this.simulateRegion(start, population);
    }

    this.props.setRegionData(regionResults);
  }

  simulateRegion(start, population) {
    let step_size = 0.1;
    let x = parseInt(1 / step_size, 10);
    let days = 200;
    let p = this.props.parameters;
    let seir_params = makeSeirParam();

    seir_params.a = 1.0 / p.incubation;
    seir_params.b = p.contact_rate;
    seir_params.g = 1 / p.infection;
    seir_params.E0 = p.e0;
    seir_params.I0 = start.Confirmed; //p.i0;
    seir_params.R0 = start.Recovered;
    seir_params.N = population;

    // TODO: replace by the actual logic
    let action_damping = calculateDamping(
      this.props.measures,
      start.date,
      days
    );

    seir_params.dampings = action_damping.map(
      (v, i) => new Damping(v.day, v.damping)
    );

    let data = simulate_seir(0, days, step_size, seir_params);

    // select only values of the days
    // TODO: this probably should be done differently
    Object.keys(data).forEach((key) => {
      data[key] = data[key].filter((v, i) => i % x === 0);
    });

    const startDate = start.date;

    const result = Object.values(
      Object.keys(data)
        .filter((k) => k !== 't')
        .reduce((acc, k) => {
          data[k].forEach((value, index) => {
            let date = new Date(startDate + index * 24 * 60 * 60 * 1000);
            date.setHours(0);
            date.setMilliseconds(0);
            date.setMinutes(0);
            date.setSeconds(0);

            date = date.getTime();
            if (!acc[date]) {
              acc[date] = {
                date
              };
            }
            acc[date] = Object.assign(acc[date], {
              [k]: parseInt(value)
            });
          });

          return acc;
        }, {})
    );

    return result;
  }

  render() {
    return (
      <Button
        onClick={this.simulate.bind(this)}
        size="sm"
        disabled={this.props.selected === null}
      >
        Simulate
      </Button>
    );
  }
}

const mapState = (state) => {
  let start = getSelectedData(state);
  if (start) {
    start = start.start;
  }
  return {
    start,
    childData: getSelectedChildData(state),
    populations: getPopulationsOfRegion(state, state.app.selected ? state.app.selected.id : 0),
    selected: state.app.selected,
    measures: getActiveMeasures(state.measures),
    parameters: getParameterMap(state.parameters)
  };
};

export default connect(mapState, { setData, setRegionData })(Simulation);
