import React, {Component} from 'react';
import {withTranslation} from 'react-i18next';
import {merge} from '../../../common/utils';

import * as moment from 'moment';
import * as numeral from 'numeral';

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

const lineProps = {
  type: 'monotone',
  dot: false,
  strokeWidth: '3'
};

const rkiLineProps = {
  strokeDasharray: '5 5',
  strokeWidth: '2'
};

const lines = [
  {
    dataKey: 'E',
    label: 'parameters.exposed',
    props: {
      stroke: '#ac58e5',
      ...lineProps
    }
  },
  {
    dataKey: 'R',
    label: 'parameters.recovered',
    props: {
      stroke: '#E0488B',
      ...lineProps
    }
  },
  {
    dataKey: 'I',
    label: 'parameters.infected',
    props: {
      stroke: '#9fd0cb',
      ...lineProps
    }
  },
  {
    dataKey: 'S',
    label: 'parameters.sus',
    props: {
      stroke: '#e0d33a',
      ...lineProps
    }
  },
  {
    dataKey: 'Recovered',
    label: 'rki.recovered',
    props: {
      stroke: '#7566ff',
      ...lineProps,
      ...rkiLineProps
    }
  },
  {
    dataKey: 'Confirmed',
    label: 'rki.infected',
    props: {
      stroke: '#533f82',
      ...lineProps,
      ...rkiLineProps
    }
  },
  {
    dataKey: 'Deaths',
    label: 'rki.deaths',
    props: {
      stroke: '#7a255d',
      ...lineProps,
      ...rkiLineProps
    }
  }
];

const longDateFormat = time => new Date(time).toLocaleDateString(undefined, {
  year: "numeric",
  month: "long",
  day: "2-digit"
});

const shortDateFormat = time => new Date(time).toLocaleDateString(undefined, {
  month: "short",
  day: "2-digit"
});


const numberFormat = (number) => {
  return numeral(number).format('0,0');
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

  constructor(props) {
    super(props);
    this.selectBar = this.selectBar.bind(this);
    this.state = {
      lines
    };
  }

  selectBar(event) {
    let updated = [];
    for (let i = 0; i < this.state.lines.length; i++) {
      let line = this.state.lines[i];
      if (line.dataKey !== event.value) {
        updated.push(line);
      } else {
        updated.push({
          ...line,
          inactive: line.inactive === undefined ? true : !line.inactive
        });
      }
    }
    this.setState({
      lines: updated
    });
  }

  translate(label) {
    const {t} = this.props;
    return t(this.state.lines.find((line) => line.dataKey === label).label);
  }

  prepareData() {
    const x = merge(
      JSON.parse(JSON.stringify(this.props.rki)),
      JSON.parse(JSON.stringify(this.props.seir)),
      'date'
    );

    try {
      x.sort(function (a, b) {
        return a.date - b.date;
      });

      return x;
    } catch (e) {
      return [];
    }
  }

  payload() {
    return this.state.lines.map((line) => {
      return {
        value: line.dataKey,
        type: 'line',
        id: line.dataKey,
        inactive: line.inactive || false,
        color: line.props.stroke
      };
    });
  }

  render() {
    const data = this.prepareData();
    const {t} = this.props;
    return (
      <ResponsiveContainer width="100%" height="80%">
        <LineChart
          data={data}
          margin={{top: 30, right: 30, left: 20, bottom: 5}}
        >
          <XAxis
            dataKey="date"
            tickFormatter={shortDateFormat}
          />
          <YAxis
            label={{
              value: t('population'),
              angle: 0,
              position: 'top',
              offset: 15
            }}
            tickFormatter={numberFormat}
          />
          <CartesianGrid strokeDasharray="3 3"/>
          <Tooltip
            offset={20}
            labelFormatter={longDateFormat}
            formatter={(value, name, index) => [
              numberFormat(value),
              this.translate(name)
            ]}
            allowEscapeViewBox={{x: true, y: false}}
            active={true}
            contentStyle={{
              backgroundColor: 'rgba(255, 255, 255, 0.8)'
            }}
            itemStyle={{
              margin: 0,
              padding: 0
            }}
          />
          {this.state.lines.map((line) => {
            return (
              <Line
                key={line.dataKey}
                dataKey={line.dataKey + (line.inactive ? ' ' : '')}
                {...line.props}
              />
            );
          })}
          <Brush
            dataKey="date"
            tickFormatter={shortDateFormat}
          />
          <Legend
            formatter={this.translate.bind(this)}
            onClick={this.selectBar}
            payload={this.payload()}
            layout="vertical"
            align="right"
            verticalAlign="top"
            wrapperStyle={{
              'padding-left': '1em'
            }}
          />
        </LineChart>
      </ResponsiveContainer>
    );
  }
}

const TranslatedChart = withTranslation()(SEIRChart);

export {TranslatedChart as SEIRChart};
