import React, { PureComponent } from 'react';
import { connect } from 'react-redux';
import { getActiveMeasures } from '../redux/measures';
import { getParameterMap } from '../redux/parameters';
import { Button } from 'reactstrap';

import { simulate_seir, makeSeirParam, Damping } from '../common/seir.js';
import { setData } from '../redux/seir';

class Simulation extends PureComponent {
  simulate() {
    let step_size = 0.1;
    let x = parseInt(1 / step_size, 10);
    let days = 200;
    let p = this.props.parameters;
    let seir_params = makeSeirParam();

    seir_params.a = 1.0 / p.incubation;
    seir_params.b = p.contact_rate;
    seir_params.g = 1 / p.infection;
    seir_params.E0 = p.e0;
    seir_params.I0 = 1000; //p.i0;
    seir_params.R0 = p.r0;
    seir_params.N = 100000;

    // TODO: replace by the actual logic
    let action_damping = [{ day: 0, damping: 1 }];

    seir_params.dampings = action_damping.map(
      (v, i) => new Damping(v.day, v.damping)
    );

    let data = simulate_seir(0, days, step_size, seir_params);

    // select only values of the days
    // TODO: this probably should be done differently
    Object.keys(data).forEach((key) => {
      data[key] = data[key].filter((v, i) => i % x === 0);
    });

    const startDate = Date.parse('2020-02-24');

    Object.keys(data)
      .filter((k) => k !== 'time')
      .forEach((k) => {
        data[k] = data[k].map((value, index) => {
          return {
            value: parseInt(value, 10),
            date: startDate + index * 24 * 60 * 60 * 1000
          };
        });
      });

    this.props.setData(data);
  }

  render() {
    return <Button onClick={this.simulate.bind(this)}>Simulate</Button>;
  }
}

const mapState = (state) => {
  return {
    measures: getActiveMeasures(state.measures),
    parameters: getParameterMap(state.parameters)
  };
};

export default connect(mapState, { setData })(Simulation);
