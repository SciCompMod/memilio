import React, { Component } from 'react';
import { connect } from 'react-redux';
import { withTranslation } from 'react-i18next';
import { Graphs as D3Graphs } from '../../common/graphs';
import { simulate_seir, makeSeirParam, Damping } from '../../common/seir.js';
import { getActiveMeasures } from '../../redux/measures';
import { getParameterMap } from '../../redux/parameters';
import * as d3 from 'd3';

import './Graphs.scss';

class OldGraphs extends Component {
  constructor(props) {
    super(props);
    this.state = {};
    this.node = React.createRef();
  }

  componentDidMount() {
    this.graphs = new D3Graphs(this.node.current);
  }

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
    console.log(data);
    // select only values of the days
    Object.keys(data).forEach((key) => {
      data[key] = data[key].filter((v, i) => i % x === 0);
    });

    let startDate = d3.timeParse('%d.%m.%Y')('24.02.2020');
    let result = [];
    // restruct data
    for (let i = 0; i < days; i++) {
      result.push({
        day: new Date(startDate.getTime() + i * 24 * 60 * 60 * 1000),
        cases: [data.S[i], data.E[i], data.I[i], data.R[i]]
      });
    }

    Object.keys(data)
      .filter((k) => k !== 'time')
      .forEach((k) => {
        data[k] = data[k].map((value, index) => {
          return {
            value,
            date: new Date(startDate.getTime() + index * 24 * 60 * 60 * 1000)
          };
        });
      });

    console.log(data);

    this.graphs.visualize(result, this.props.measures);
  }

  render() {
    if (this.graphs) {
      this.simulate();
    }
    return (
      <div className="graphs">
        <svg ref={this.node}></svg>
      </div>
    );
  }
}

const mapState = (state) => {
  return {
    measures: getActiveMeasures(state.measures),
    parameters: getParameterMap(state.parameters)
  };
};

export default connect(mapState, null)(withTranslation()(OldGraphs));
