import React, {Component} from 'react';
import {withTranslation} from 'react-i18next';
import * as am4core from '@amcharts/amcharts4/core';
import * as am4charts from '@amcharts/amcharts4/charts';

import './InfectionChart.scss';

/**
 * This component wraps an amChart to display infection related data of RKI and SEIR simulations.
 */
class InfectionChart extends Component {
  static defaultProps = {
    seir: {
      S: [],
      E: [],
      I: [],
      R: [],
    },
    rki: [],
    measures: [],
  };

  /** @type XYChart */
  #chart = null;

  /** @type ValueAxis */
  #yAxis = null;

  /** @type Map<string, LineSeries> */
  #series = new Map();

  componentDidMount() {
    this.#chart = am4core.create('dataChartDiv', am4charts.XYChart);
    this.#chart.xAxes.push(new am4charts.DateAxis());
    this.#yAxis = this.#chart.yAxes.push(new am4charts.ValueAxis());

    this.#chart.legend = new am4charts.Legend();
    this.#chart.cursor = new am4charts.XYCursor();
    this.#chart.cursor.maxTooltipDistance = -1;

    this.#chart.dateFormatter.dateFormat = {
      month: 'short',
      day: '2-digit',
    };

    this.#chart.scrollbarX = new am4charts.XYChartScrollbar();
    this.#chart.scrollbarX.parent = this.#chart.bottomAxesContainer;
    this.#chart.data = [];

    const {t} = this.props;
    this.createSeries(t('rki.infected'), 'rkic');
    this.createSeries(t('rki.deaths'), 'rkid');
    this.createSeries(t('rki.recovered'), 'rkir');

    this.createSeries(t('parameters.sus'), 'seirs');
    this.createSeries(t('parameters.exposed'), 'seire');
    this.createSeries(t('parameters.infected'), 'seiri');
    this.createSeries(t('parameters.recovered'), 'seirr');

    this.updateData();
  }

  /**
   * Creates a new series and adds it to the chart.
   * @private
   * @param name{string} The name to be displayed.
   * @param valueName{string} The name of the value field.
   */
  createSeries(name, valueName) {
    const series = this.#chart.series.push(new am4charts.LineSeries());
    series.name = name;
    series.dataFields.valueY = valueName;
    series.dataFields.dateX = 'date';

    series.tooltipText = `${name}: {valueY}`;
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

    this.#chart.scrollbarX.series.push(series);
    this.#series.set(valueName, series);
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (prevProps.rki !== this.props.rki || prevProps.seir !== this.props.seir) {
      this.updateData();
    }

    if (prevProps.logChart !== this.props.logChart) {
      this.#yAxis.logarithmic = this.props.logChart;

      // logarithmic scale can't handle 0 values, so we set them to 1, which will be displayed as 0.
      this.#yAxis.treatZeroAs = this.props.logChart ? 1 : 0;
    }
  }

  /**
   * Transforms the incoming data into a format that the chart needs and updates
   * the chart and series accordingly.
   * @private
   */
  updateData() {
    const timeMap = new Map();

    if (this.props.rki === null || Object.keys(this.props.rki).length === 0) {
      this.#series.get('rkic').hide();
      this.#series.get('rkid').hide();
      this.#series.get('rkir').hide();
    } else {
      this.#series.get('rkic').show();
      this.#series.get('rkid').show();
      this.#series.get('rkir').show();
    }

    if (this.props.seir === null || Object.keys(this.props.seir).length === 0) {
      this.#series.get('seirs').hide();
      this.#series.get('seire').hide();
      this.#series.get('seiri').hide();
      this.#series.get('seirr').hide();
    } else {
      this.#series.get('seirs').show();
      this.#series.get('seire').show();
      this.#series.get('seiri').show();
      this.#series.get('seirr').show();
    }

    if (this.props.rki !== null) {
      for (let entry of Object.values(this.props.rki)) {
        if (!timeMap.has(entry.date)) {
          timeMap.set(entry.date, {});
        }

        const curr = timeMap.get(entry.date);
        curr.rkic = entry.Confirmed;
        curr.rkid = entry.Deaths;
        curr.rkir = entry.Recovered;
      }
    }

    if (this.props.seir !== null) {
      for (let entry of Object.values(this.props.seir)) {
        if (!timeMap.has(entry.date)) {
          timeMap.set(entry.date, {});
        }

        const curr = timeMap.get(entry.date);
        curr.seirs = entry.S;
        curr.seire = entry.E;
        curr.seiri = entry.I;
        curr.seirr = entry.R;
      }
    }

    const data = [];
    for (let [key, value] of timeMap.entries()) {
      data.push({
        date: key,
        ...value,
      });
    }

    data.sort((a, b) => a.date - b.date);
    this.#chart.data = data;
  }

  render() {
    return (
      <div style={{height: '100%'}}>
        <div id="dataChartDiv" style={{height: '100%'}} />
      </div>
    );
  }
}

const TranslatedChart = withTranslation()(InfectionChart);

export {TranslatedChart as InfectionChart};
