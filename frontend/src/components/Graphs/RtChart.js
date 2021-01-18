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
    dataset: '',
  };

  /** @type XYChart */
  #chart = null;

  /** @type DateAxis */
  #xAxis = null;

  /** @type ValueAxis */
  #yAxis = null;

  #legend = null;

  /** @type Map<string, LineSeries> */
  #series = null;

  componentDidMount() {
    this.init();
  }

  init() {
    this.#chart = am4core.create(`rt-chart-${this.props.id}`, am4charts.XYChart);
    this.#xAxis = this.#chart.xAxes.push(new am4charts.DateAxis());
    this.#yAxis = this.#chart.yAxes.push(new am4charts.ValueAxis());

    // disable axis tooltips
    this.#xAxis.cursorTooltipEnabled = false;
    this.#yAxis.cursorTooltipEnabled = false;

    // add legend
    this.#legend = new am4charts.Legend();
    this.#chart.legend = this.#legend;

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
    this.#series = new Map();
    for (let i = 0; i < this.props.series.length; i++) {
      this.createSeries(this.props.series[i]);
    }

    // set data
    this.#chart.data = this.props.data;
  }

  componentDidUpdate(prevProps) {
    if (prevProps.id !== this.props.id) {
      this.#chart.dispose();
      this.init();
    }

    // show / hide series
    for (let i = 0; i < this.props.series.length; i++) {
      const s = this.props.series[i];
      let serie = this.#series.get(s.key);
      if (!serie) {
        return;
      }
      console.log(s, serie);
      console.log(s.key, s.isHidden);
      if (s.isHidden) {
        // check if already hidden and hide if not

        // if series is still visible remove it
        console.log('hide series', s.key);
        //this.#chart.series.removeIndex(index);
        //this.#series.delete(s.key);
        serie.hide();
      } else {
        // check if already visible and show if not

        // if series is not visible add it
        console.log('show series', s.key);
        serie.show();
        //this.createSeries(s);
      }
    }

    // we only need to update the value when the district value has changes
    if (this.props.district !== prevProps.district) {
      this.#chart.data = this.props.data;
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
    const series = new am4charts.LineSeries();

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

    series.autoDispose = false;

    this.#chart.series.push(series);
    if (isHidden) {
      series.hide();
    }

    this.#series.set(key, series);
  }

  render() {
    return <div id={`rt-chart-${this.props.id}`} style={{width: '100%', height: '100%'}}></div>;
  }
}
