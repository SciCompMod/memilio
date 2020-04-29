import React, { Component, PureComponent } from 'react';
import { connect } from 'react-redux';
import { withTranslation } from 'react-i18next';
import { scaleTime } from 'd3-scale';
import { timeFormat } from 'd3-time-format';

import ResponsiveXYFrame from 'semiotic/lib/ResponsiveXYFrame';
import XYFrame from 'semiotic/lib/XYFrame';

import './SEIRChart.scss';

const translationMap = {
  E: 'parameters.exposed',
  R: 'parameters.recovered',
  I: 'parameters.infected',
  S: 'parameters.sus'
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

const tooltipStyles = {
  header: {
    fontWeight: 'bold',
    borderBottom: 'thin solid black',
    marginBottom: '10px',
    textAlign: 'center'
  },
  lineItem: { position: 'relative', display: 'block', textAlign: 'left' },
  title: { display: 'inline-block', margin: '0 5px 0 15px' },
  value: { display: 'inline-block', fontWeight: 'bold', margin: '0' },
  wrapper: {
    background: 'rgba(255,255,255,0.8)',
    minWidth: 'max-content',
    whiteSpace: 'nowrap'
  }
};

const plotProps = {
  size: [500, 300],
  margin: { left: 80, bottom: 90, right: 10, top: 100 },
  xScaleType: scaleTime(),
  xAccessor: (d) => new Date(d.date),
  yAccessor: (d) => d.value,
  yExtent: [0],
  responsiveWidth: true,
  pointStyle: { fill: 'none', stroke: 'black', strokeWidth: '1.5px' },
  className: 'sharedTooltip',
  axes: [
    {
      orient: 'left',
      label: 'Personen',
      tickFormat: function (e) {
        if (e > 1000) {
          return e / 1000 + 'k';
        }
        return e;
      }
    },
    {
      orient: 'bottom',
      ticks: 6,
      tickFormat: (d) => timeFormat('%d.%m')(d),
      label: { name: 'Zeit', locationDistance: 55 }
    }
  ],
  lineStyle: (d, i) => ({
    stroke: theme[i],
    strokeWidth: 2,
    fill: 'none'
  }),
  hoverAnnotation: [
    { type: 'x', disable: ['connector', 'note'] },
    { type: 'frame-hover' }
  ],
  annotations: [],
  annotationSettings: {
    layout: {
      type: 'marginalia',
      orient: 'top',
      marginOffset: 20,
      padding: 10,
      lineHeight: 15,
      characterWidth: 12
    }
  }
};

class SEIRChart extends Component {
  static defaultProps = {
    S: [],
    E: [],
    I: [],
    R: [],
    measures: []
  };

  state = { extend: null };

  prepareProps() {
    const { t } = this.props;
    const lines = ['S', 'E', 'I', 'R'].map((k) => {
      return {
        title: t(translationMap[k]),
        coordinates: this.props[k]
      };
    });

    const annotations = this.props.measures
      .map((m, i) => {
        return m.intervals.map((interval) => ({
          className: `measure measure-${i}`,
          type: 'bounds',
          bounds: [
            { date: new Date(interval.start) },
            { date: new Date(interval.end) }
          ],
          label: t(m.label)
        }));
      })
      .flat();

    return {
      ...plotProps,
      lines,
      annotations,
      tooltipContent: (d) => {
        const points = lines.map((line, index) => {
          return {
            color: theme[index],
            title: line.title,
            ...line.coordinates.find((i) => {
              // Search the lines for a similar x value for vertical shared tooltip
              // Can implement a 'close enough' conditional here too (fuzzy equality)
              return i.date === d.x.getTime();
            })
          };
        });

        const returnArray = [
          <div key={'header_multi'} style={tooltipStyles.header}>
            {`${timeFormat('%d.%m.%Y')(new Date(d.x))}`}
          </div>
        ];

        points.forEach((point, i) => {
          const title = point.title;
          const valString = `${point.value}`;

          returnArray.push([
            <div key={`tooltip_line_${i}`} style={tooltipStyles.lineItem}>
              <p
                key={`tooltip_color_${i}`}
                style={{
                  width: '10px',
                  height: '10px',
                  backgroundColor: point.color,
                  display: 'inline-block',
                  position: 'absolute',
                  top: '8px',
                  left: '0',
                  margin: '0'
                }}
              />
              <p
                key={`tooltip_p_${i}`}
                style={tooltipStyles.title}
              >{`${title}:`}</p>
              <p key={`tooltip_p_val_${i}`} style={tooltipStyles.value}>
                {valString}
              </p>
            </div>
          ]);
        });

        return (
          <div className="tooltip-content" style={tooltipStyles.wrapper}>
            {returnArray}
          </div>
        );
      }
    };
  }

  renderFrame() {
    if (this.props.E.length > 0) {
      const props = this.prepareProps();
      const { annotations, annotationSettings, ...brushProps } = props;

      brushProps.size = [500, 200];
      brushProps.margin = { ...brushProps.margin, top: 0, bottom: 10 };

      return (
        <>
          <div className="plot">
            <ResponsiveXYFrame {...props} />
          </div>
          <div className="brush">
            <ResponsiveXYFrame
              {...brushProps}
              interaction={{
                end: (e) => {
                  console.log(e);
                  this.setState({ extent: e });
                },
                brush: 'xBrush',
                extent: this.state.extent
              }}
            />
          </div>
        </>
      );
    }

    return <></>;
  }

  render() {
    return this.renderFrame();
  }
}

export default withTranslation()(SEIRChart);
