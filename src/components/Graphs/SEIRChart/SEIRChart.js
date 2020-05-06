import React, { Component } from 'react';
import { withTranslation } from 'react-i18next';
import { merge } from '../../../common/utils';

import * as moment from 'moment';

import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  Brush,
  ResponsiveContainer
} from 'recharts';

import './SEIRChart.scss';

const translationMap = {
  E: 'parameters.exposed',
  R: 'parameters.recovered',
  I: 'parameters.infected',
  S: 'parameters.sus',
  AnzahlGenesen: 'rki.recovered',
  AnzahlFall: 'rki.infected',
  AnzahlTodesfall: 'rki.deaths'
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
    seir: {
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
    const x = merge(
      JSON.parse(JSON.stringify(this.props.rki)),
      JSON.parse(JSON.stringify(this.props.seir)),
      'date'
    );
    x.sort(function (a, b) {
      return a.date - b.date;
    });

    return x;
  }

  render() {
    const data = this.prepareData();
    return (
      <ResponsiveContainer width="100%" height="80%">
        <LineChart
          data={data}
          margin={{ top: 5, right: 30, left: 20, bottom: 5 }}
        >
          <XAxis dataKey="date" tickFormatter={dateFormat} />
          <YAxis />
          <CartesianGrid strokeDasharray="3 3" />
          <Tooltip
            labelFormatter={dateFormat}
            formatter={(value, name, index) => [value, this.translate(name)]}
          />

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
          <Line dataKey="S" stroke={theme[4]} {...lineProps} />
          <Line dataKey="E" stroke={theme[3]} {...lineProps} />
          <Line dataKey="I" stroke={theme[0]} {...lineProps} />
          <Line dataKey="R" stroke={theme[2]} {...lineProps} />
          <Brush dataKey="date" tickFormatter={dateFormat} />
          <Legend formatter={this.translate.bind(this)} />
        </LineChart>
      </ResponsiveContainer>
    );
  }
}

const TranslatedChart = withTranslation()(SEIRChart);

export { TranslatedChart as SEIRChart };
