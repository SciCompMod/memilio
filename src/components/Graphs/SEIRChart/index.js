import React, { Component, PureComponent } from 'react';
import { connect } from 'react-redux';
import { withTranslation } from 'react-i18next';
import { scaleTime } from 'd3-scale';
import { timeFormat } from 'd3-time-format';
import * as d3_time from 'd3-time';
import * as d3_scale from 'd3-scale';

import * as moment from 'moment';

import ResponsiveXYFrame from 'semiotic/lib/ResponsiveXYFrame';
import XYFrame from 'semiotic/lib/XYFrame';

import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  Brush
} from 'recharts';

import './SEIRChart.scss';

const keyMap = {
  date: 'Meldedatum',
  I: 'AnzahlFall',
  D: 'AnzahlTodesfall',
  R: 'AnzahlGenesen'
};

const translationMap = {
  E: 'parameters.exposed',
  R: 'parameters.recovered',
  I: 'parameters.infected',
  S: 'parameters.sus',
  E: 'parameters.exposed',
  AnzahlGenesen: 'parameters.recovered',
  AnzahlFall: 'parameters.infected',
  AnzahlTodesfall: 'parameters.deaths'
};

const theme = [
  '#ac58e5',
  '#E0488B',
  '#9fd0cb',
  '#e0d33a',
  '#7566ff',
  '#533f82',
  '#7a255d',
  '#365350',
  '#a19a11',
  '#3f4482'
];

const lineProps = {
  type: 'monotone',
  dot: false,
  strokeWidth: '3'
};

const rkiLineProps = {
  strokeDasharray: '5 5',
  strokeWidth: '2'
};

const dateFormat = (time) => {
  return moment(time).format('DD. MMM');
};

class SEIRChart extends Component {
  static defaultProps = {
    simulation: {
      S: [],
      E: [],
      I: [],
      R: []
    },
    rki: [],
    measures: []
  };

  translate(label) {
    const { t } = this.props;
    return t(translationMap[label]);
  }

  prepareData() {
    return this.props.rki;
  }

  render() {
    const { t } = this.props;
    const data = this.prepareData();
    return (
      <LineChart
        width={600}
        height={300}
        data={data}
        margin={{ top: 5, right: 30, left: 20, bottom: 5 }}
      >
        <XAxis dataKey="Meldedatum" tickFormatter={dateFormat} />
        <YAxis />
        <CartesianGrid strokeDasharray="3 3" />
        <Tooltip
          labelFormatter={dateFormat}
          formatter={(value, name, index) => [value, this.translate(name)]}
        />
        <Legend formatter={this.translate.bind(this)} />
        <Line
          dataKey="AnzahlFall"
          stroke={theme[0]}
          {...lineProps}
          {...rkiLineProps}
        />
        <Line
          dataKey="AnzahlTodesfall"
          stroke={theme[1]}
          {...lineProps}
          {...rkiLineProps}
        />
        <Line
          dataKey="AnzahlGenesen"
          stroke={theme[2]}
          {...lineProps}
          {...rkiLineProps}
        />
        <Brush data={data} dataKey="Meldedatum" tickFormatter={dateFormat} />
      </LineChart>
    );
  }
}

export default withTranslation()(SEIRChart);
