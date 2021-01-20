import * as am4core from '@amcharts/amcharts4/core';
import * as am4charts from '@amcharts/amcharts4/charts';
import * as am4maps from '@amcharts/amcharts4/maps';
import {Subject} from 'rxjs';

export default class InteractiveHeatMap extends Subject {
  /** @type MapChart */
  #map = null;

  /** @type MapPolygonSeries */
  #countySeries = null;

  /** @type HeatLegend */
  #countyHeatLegend = null;

  #min = null;
  #max = null;

  isReady = false;

  constructor(id, options) {
    super();
    const {showExport, showLegend} = Object.assign({showExport: false, showLegend: true}, options);

    this.disable_map_hit = false;

    this.#map = am4core.create(id, am4maps.MapChart);

    this.#map.events.on('ready', () => {
      this.next('ready');
      this.isReady = true;
    });

    this.#map.projection = new am4maps.projections.Mercator();

    if (showExport) {
      this.#map.exporting.menu = new am4core.ExportMenu();
      this.#map.exporting.menu.align = 'right';
      this.#map.exporting.menu.verticalAlign = 'bottom';
      this.#map.hiddenState.transitionDuration = 100;
      this.#map.maxPanOut = 0.0;
    }

    this.#map.numberFormatter.numberFormat = '#,###.##';

    this.#countySeries = new am4maps.MapPolygonSeries();

    const heatColors = [
      am4core.color('rgb(255, 255, 255)'),
      am4core.color('rgb(51, 204, 51)'),
      am4core.color('rgb(230, 204, 51)'),
      am4core.color('rgb(230, 153, 0)'),
      am4core.color('rgb(230, 0, 0)'),
      am4core.color('rgb(204, 0, 204)'),
    ];

    const countyPolygonTemplate = this.#countySeries.mapPolygons.template;
    countyPolygonTemplate.tooltipText = '{BEZ} {GEN}: {value}';
    countyPolygonTemplate.nonScalingStroke = true;
    countyPolygonTemplate.applyOnClones = true;
    countyPolygonTemplate.stroke = am4core.color('rgb(221, 221, 221)');
    countyPolygonTemplate.fill = am4core.color('rgb(255, 255, 255)');

    if (showLegend) {
      this.#countyHeatLegend = this.#map.createChild(am4charts.HeatLegend);
      this.#countyHeatLegend.series = this.#countySeries;
      this.#countyHeatLegend.height = am4core.percent(20);
      this.#countyHeatLegend.align = 'left';
      this.#countyHeatLegend.valign = 'bottom';
      this.#countyHeatLegend.orientation = 'vertical';
      this.#countyHeatLegend.marginBottom = 50;
      this.#countyHeatLegend.minValue = 0;
      this.#countyHeatLegend.maxValue = 1;
      this.#countyHeatLegend.fontSize = '.8rem';

      // Allow the heatLegend to function in general
      this.#countyHeatLegend.minColor = heatColors[0];
      this.#countyHeatLegend.maxColor = heatColors[heatColors.length - 1];

      // This event ensures, that the legend works, when the window resizes.
      // It looks like it does nothing, but it actually works -.-
      this.#countyHeatLegend.markers.template.events.on('sizechanged', (event) => {
        // eslint-disable-next-line no-self-assign
        return (event.target.fill = event.target.fill);
      });

      // Override heatLegend gradient
      const gradient = new am4core.LinearGradient();
      gradient.rotation = -90;
      heatColors.forEach((color) => gradient.addColor(color));

      this.#countyHeatLegend.markers.template.adapter.add('fill', () => gradient);
    }
    countyPolygonTemplate.events.on('hit', (e) => {
      const item = e.target.dataItem.dataContext;
      this.next({rs: item.RS, gen: item.GEN, bez: item.BEZ});
      this.disable_map_hit = true;
    });

    this.#map.events.on('hit', (e) => {
      if (this.disable_map_hit) {
        this.disable_map_hit = false;
        return;
      }
      this.next('reset');
    });

    const series = this.#map.series.push(this.#countySeries);

    const getColor = (workingValue, minValue, maxValue) => {
      // calculate percentage and restrict it between 0 and 1
      const percent = Math.max(0, Math.min(1, (workingValue - minValue) / (maxValue - minValue)));
      const intervals = heatColors.length - 1;
      const fract = 1 / intervals;

      if (am4core.type.isNumber(percent)) {
        const colorIndex = Math.max(0, parseInt(Math.ceil(intervals * percent - 1), 10));

        return new am4core.Color(
          am4core.colors.interpolate(
            heatColors[colorIndex].rgb,
            heatColors[colorIndex + 1].rgb,
            (percent - colorIndex * fract) / fract
          )
        );
      }

      return am4core.color('white');
    };

    // for some reason this event is triggered multiple times,
    // therfore we only apply the new colors only when the max value changed!
    series.events.on('validated', (event) => {
      if (event.target.children.values.length > 2) {
        const {high} = this.#countySeries.dataItem.values['value'];

        const min = this.#min || 0;
        const max = this.#max || high;

        if (showLegend) {
          this.#countyHeatLegend.minValue = min;
          this.#countyHeatLegend.maxValue = max;
        }

        event.target.children.values.forEach((mapPolygon, index) => {
          if (index > 1) {
            const {workingValue} = mapPolygon.dataItem.values['value'];
            const color = getColor(workingValue, min, max);
            mapPolygon.fill = color;
          }
        });
      }
    });

    this.#countySeries.useGeodata = true;
    this.#countySeries.geodataSource.url = 'assets/counties/lk_germany_reduced.geojson';

    this.#map.series.push(this.#countySeries);
    //this.#map.invalidate();
  }

  /** @param values {Map<string, number>} */
  async setValues(values) {
    if (!values) {
      return;
    }

    for (let countyDatum of this.#countySeries.data) {
      countyDatum.value = values.get(countyDatum.RS) ?? 0;
    }

    this.#countySeries.invalidateRawData();
  }

  setLegendMinMax(min, max) {
    // ensure min max are numbers
    min = parseFloat(min);
    max = parseFloat(max);

    if (isNaN(min)) {
      min = 0;
    }

    if (isNaN(max)) {
      max = undefined;
    }
    this.#min = min;
    this.#max = max;

    this.#countySeries.invalidateRawData();
  }
}
