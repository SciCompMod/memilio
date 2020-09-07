import React, {Component} from 'react';
import {withTranslation} from 'react-i18next';
import * as am4core from '@amcharts/amcharts4/core';
import * as am4charts from '@amcharts/amcharts4/charts';

import './Chart.scss';

/**
 * This component wraps an amChart to display infection related data of RKI and SEIR simulations.
 */
class Chart extends Component {
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

  /** @private
   *  @type XYChart */
  _chart = null;

  /** @private
   *  @type ValueAxis */
  _yAxis = null;

  /** @private
   *  @type Map<string, LineSeries> */
  _series = new Map();

  componentDidMount() {
    this._chart = am4core.create('dataChartDiv', am4charts.XYChart);
    this._chart.xAxes.push(new am4charts.DateAxis());
    this._yAxis = this._chart.yAxes.push(new am4charts.ValueAxis());

    this._chart.legend = new am4charts.Legend();
    this._chart.cursor = new am4charts.XYCursor();
    this._chart.cursor.maxTooltipDistance = -1;

    this._chart.dateFormatter.dateFormat = {
      month: 'short',
      day: '2-digit',
    };

    this._chart.scrollbarX = new am4charts.XYChartScrollbar();
    this._chart.scrollbarX.parent = this._chart.bottomAxesContainer;
    this._chart.data = [];

    const {t} = this.props;
    this._createSeries(t('rki.infected'), 'rkic');
    this._createSeries(t('rki.deaths'), 'rkid');
    this._createSeries(t('rki.recovered'), 'rkir');

    this._createSeries(t('parameters.sus'), 'seirs');
    this._createSeries(t('parameters.exposed'), 'seire');
    this._createSeries(t('parameters.infected'), 'seiri');
    this._createSeries(t('parameters.recovered'), 'seirr');

    this._updateData();
  }

  /**
   * Creates a new series and adds it to the chart.
   * @private
   * @param name{string} The name to be displayed.
   * @param valueName{string} The name of the value field.
   */
  _createSeries(name, valueName) {
    const series = this._chart.series.push(new am4charts.LineSeries());
    series.name = name;
    series.dataFields.valueY = valueName;
    series.dataFields.dateX = 'date';

    series.tooltipText = `${name}: {valueY}`;
    series.tooltip.getFillFromObject = false;
    series.tooltip.background.fill = am4core.color('#fff');
    series.tooltip.label.fill = am4core.color('#00');
    series.adapter.add('tooltipText', (e) => {
      let text = '[bold]{dateX}[]\n';
      this._chart.series.each((item) => {
        if (item.visible) {
          text += `[${item.stroke.hex}]●[/] ${item.name}: {${item.dataFields.valueY}}\n`;
        }
      });
      return text;
    });

    this._chart.scrollbarX.series.push(series);
    this._series.set(valueName, series);
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (prevProps.rki !== this.props.rki || prevProps.seir !== this.props.seir) {
      this._updateData();
    }

    if (prevProps.logChart !== this.props.logChart) {
      this._yAxis.logarithmic = this.props.logChart;

      // logarithmic scale can't handle 0 values, so we set them to 1, which will be displayed as 0.
      this._yAxis.treatZeroAs = this.props.logChart ? 1 : 0;
    }
  }

  /**
   * Transforms the incoming data into a format that the chart needs and updates
   * the chart and series accordingly.
   * @private
   */
  _updateData() {
    const timeMap = new Map();

    if (this.props.rki === null || Object.keys(this.props.rki).length === 0) {
      this._series.get('rkic').hide();
      this._series.get('rkid').hide();
      this._series.get('rkir').hide();
    } else {
      this._series.get('rkic').show();
      this._series.get('rkid').show();
      this._series.get('rkir').show();
    }

    if (this.props.seir === null || Object.keys(this.props.seir).length === 0) {
      this._series.get('seirs').hide();
      this._series.get('seire').hide();
      this._series.get('seiri').hide();
      this._series.get('seirr').hide();
    } else {
      this._series.get('seirs').show();
      this._series.get('seire').show();
      this._series.get('seiri').show();
      this._series.get('seirr').show();
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
    this._chart.data = data;
  }

  render() {
    return (
      <div style={{height: '100%'}}>
        <div id="dataChartDiv" style={{height: '100%'}} />
      </div>
    );
  }
}

const TranslatedChart = withTranslation()(Chart);

export {TranslatedChart as Chart};
