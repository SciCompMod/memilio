import React, {Component} from 'react';

import * as am4core from '@amcharts/amcharts4/core';
import * as am4charts from '@amcharts/amcharts4/charts';

/**
 * This chart displays the reproduction values RT and RT relative
 */
export default class RtChart extends Component {
  static defaultProps = {
    data: {},
    district: '',
  };

  /** @type XYChart */
  #chart = null;

  /** @type DateAxis */
  #xAxis = null;

  /** @type ValueAxis */
  #yAxis = null;

  /** @type Map<string, LineSeries> */
  #series = new Map();

  componentDidMount() {
    this.#chart = am4core.create('rt-chart', am4charts.XYChart);
    this.#xAxis = this.#chart.xAxes.push(new am4charts.DateAxis());
    this.#yAxis = this.#chart.yAxes.push(new am4charts.ValueAxis());

    // disable axis tooltips
    this.#xAxis.cursorTooltipEnabled = false;
    this.#yAxis.cursorTooltipEnabled = false;

    // add xy cursor
    this.#chart.cursor = new am4charts.XYCursor();
    this.#chart.cursor.maxTooltipDistance = -1;
    this.#chart.cursor.lineY.disabled = true;

    // set number format
    this.#chart.numberFormatter.numberFormat = '#,###.##';

    // define data formats
    this.#chart.dateFormatter.dateFormat = 'dd. MMMM';
    this.#xAxis.dateFormats.setKey('day', 'dd. MMMM');
    this.#xAxis.dateFormats.setKey('month', 'MMMM');

    // create given series
    for (let i = 0; i < this.props.series.length; i++) {
      this.createSeries(this.props.series[i]);
    }

    // set data
    this.#chart.data = this.props.data;
  }

  componentDidUpdate(prevProps) {
    // we only need to update the value when the district value has changes
    if (this.props.district !== prevProps.district) {
      this.#chart.data = this.props.data;

      // show / hide series
      this.props.series.forEach((serie) => {
        if (serie.isHidden) {
          this.#series.get(serie.key).hide();
        } else {
          this.#series.get(serie.key).show();
        }
      });
    }
  }

  componentWillUnmount() {
    if (this.#chart) {
      this.#chart.dispose();
    }
  }

  /**
   * Creates a new series and adds it to the chart.
   * @private
   * @param name{string} The name to be displayed.
   * @param valueName{string} The name of the value field.
   */
  createSeries(options) {
    const {key, label, isHidden} = options;
    const series = this.#chart.series.push(new am4charts.LineSeries());
    series.name = label;
    series.dataFields.valueY = key;
    series.dataFields.dateX = 'timestamp';

    series.tooltipText = `${label}: {valueY}`;
    series.tooltip.getFillFromObject = false;
    series.tooltip.background.fill = am4core.color('#fff');
    series.tooltip.label.fill = am4core.color('#00');
    series.adapter.add('tooltipText', (e) => {
      let text = '[bold]{dateX}[]\n';
      this.#chart.series.each((item) => {
        if (item.visible) {
          text += `[${item.stroke.hex}]●[/] ${item.name}: {${item.dataFields.valueY}}\n`;
        }
      });
      return text;
    });
    series.tensionX = 0.8;
    series.tensionY = 0.8;
    series.strokeWidth = 3;

    series.hidden = isHidden || false;

    this.#series.set(key, series);
  }

  render() {
    return <div id="rt-chart" style={{width: '100%', height: '100%'}}></div>;
  }
}
