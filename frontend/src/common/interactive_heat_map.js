import * as am4core from '@amcharts/amcharts4/core';
import * as am4charts from '@amcharts/amcharts4/charts';
import * as am4maps from '@amcharts/amcharts4/maps';

export default class InteractiveHeatMap {
  /** @type MapChart */
  #map = null;

  /** @type MapPolygonSeries */
  #stateChoroplethSeries = null;

  /** @type MapPolygonSeries */
  #stateBackgroundSeries = null;

  /** @type Map<number, MapPolygonSeries> */
  #countySeries = new Map();

  /** @type HeatLegend */
  #stateHeatLegend = null;

  /** @type Map<number, HeatLegend> */
  #countyHeatLegends = new Map();

  /** @type am4core.Label */
  #dataSetLabel = null;

  /** @type IHeatRule */
  #stateHeatRule = null;

  /** @type boolean */
  #seriesHit = false;

  /** @type number */
  selectedState = -1;

  /** @type number */
  selectedCounty = -1;

  /** @callback StateSelectedCallback
   *  @param newState {StateMetaData} */
  onStateSelected = (newState) => {};

  /** @callback CountySelectedCallback
   *  @param newCounty {any} */
  onCountySelected = (newCounty) => {};

  constructor(id) {
    this.#map = am4core.create(id, am4maps.MapChart);
    this.#map.projection = new am4maps.projections.Mercator();
    this.#map.hiddenState.transitionDuration = 100;
    this.#map.maxPanOut = 0.0;

    this.#dataSetLabel = this.#map.createChild(am4core.Label);
    this.#dataSetLabel.align = 'right';
    this.#dataSetLabel.valign = 'bottom';
    this.#dataSetLabel.marginBottom = 20;
    this.#dataSetLabel.marginRight = 20;
    this.#dataSetLabel.text = 'Dataset: RKI';

    this.#stateBackgroundSeries = new am4maps.MapPolygonSeries();
    this.#stateBackgroundSeries.geodataSource.url = 'assets/german-states.geojson';
    this.#stateBackgroundSeries.useGeodata = true;
    this.#map.series.push(this.#stateBackgroundSeries);
    this.#stateBackgroundSeries.mapPolygons.template.fill = am4core.color('#AAA');
    this.#stateBackgroundSeries.mapPolygons.template.tooltipText = '{name}';
    this.#stateBackgroundSeries.hide();

    this.#stateChoroplethSeries = new am4maps.MapPolygonSeries();
    this.#stateChoroplethSeries.geodataSource.url = 'assets/german-states.geojson';
    this.#stateChoroplethSeries.useGeodata = true;
    this.#map.series.push(this.#stateChoroplethSeries);

    this.#stateHeatRule = {
      property: 'fill',
      target: this.#stateChoroplethSeries.mapPolygons.template,
      min: am4core.color('#DDD', 1),
      max: am4core.color('#F00', 1),
      logarithmic: false,
    };
    this.#stateChoroplethSeries.heatRules.push(this.#stateHeatRule);

    this.#stateHeatLegend = this.#map.createChild(am4charts.HeatLegend);
    this.#stateHeatLegend.series = this.#stateChoroplethSeries;
    this.#stateHeatLegend.height = am4core.percent(90);
    this.#stateHeatLegend.align = 'left';
    this.#stateHeatLegend.valign = 'middle';
    this.#stateHeatLegend.orientation = 'vertical';
    this.#stateHeatLegend.marginLeft = 20;

    const statePolygonTemplate = this.#stateChoroplethSeries.mapPolygons.template;
    statePolygonTemplate.tooltipText = '{name}: {value}';
    statePolygonTemplate.nonScalingStroke = true;
    statePolygonTemplate.applyOnClones = true;

    const selectStateEvent = (e) => {
      this.#seriesHit = true;
      const item = e.target.dataItem.dataContext;
      this.onStateSelected(item);
      this.onCountySelected(null);
      this.selectedCounty = -1;
      this.stateSelected(item);
    };

    statePolygonTemplate.events.on('hit', selectStateEvent);
    this.#stateBackgroundSeries.mapPolygons.template.events.on('hit', selectStateEvent);

    const zoomOutEvent = () => {
      this.onStateSelected(null);
      this.onCountySelected(null);
      this.selectedCounty = -1;
      this.stateSelected({id: -1});
    };

    this.#map.events.on('hit', () => {
      if (!this.#seriesHit) {
        zoomOutEvent();
      }
      this.#seriesHit = false;
    });

    this.#map.events.on('zoomlevelchanged', () => {
      if (this.#map.zoomLevel === 1) {
        zoomOutEvent();
      }
    });

    for (let state of STATE_METADATA) {
      this.createCountySeries(state);
    }

    this.#map.invalidate();
  }

  setDataSetName(name) {
    this.#dataSetLabel.text = 'Dataset: ' + name;
  }

  /** @param values {Map<number, number>} */
  setStateValues(values) {
    if (values) {
      for (let stateDatum of this.#stateChoroplethSeries.data) {
        stateDatum.value = values.get(stateDatum.id) ?? 0;
      }

      this.#stateChoroplethSeries.invalidateRawData();
    }
  }

  /** @param values {Map<string, number>} */
  setCountyValues(values) {
    if (!values) {
      return;
    }

    if (this.#countySeries.has(this.selectedState)) {
      const series = this.#countySeries.get(this.selectedState);
      for (let countyDatum of series.data) {
        countyDatum.value = values.get(parseInt(countyDatum.RS, 10)) ?? 0;
      }

      series.invalidateRawData();
    }
  }

  /** @private
   *  @param state {StateMetaData} */
  createCountySeries(state) {
    const newSeries = new am4maps.MapPolygonSeries();
    newSeries.geodataSource.url = 'assets/counties/' + state.fileName + '.geojson';
    newSeries.useGeodata = true;
    newSeries.hidden = true;

    const countyPolygonTemplate = newSeries.mapPolygons.template;
    countyPolygonTemplate.tooltipText = '{id}: {value}';
    countyPolygonTemplate.nonScalingStroke = true;
    countyPolygonTemplate.applyOnClones = true;

    const countyHeatRule = {
      property: 'fill',
      target: countyPolygonTemplate,
      min: am4core.color('#DDD', 1),
      max: am4core.color('#F00', 1),
      logarithmic: false,
    };

    newSeries.heatRules.push(countyHeatRule);

    const heatLegend = this.#map.createChild(am4charts.HeatLegend);
    heatLegend.series = newSeries;
    heatLegend.height = am4core.percent(90);
    heatLegend.align = 'left';
    heatLegend.valign = 'middle';
    heatLegend.orientation = 'vertical';
    heatLegend.marginLeft = 20;

    countyPolygonTemplate.events.on('hit', (e) => {
      this.#seriesHit = true;
      const item = e.target.dataItem.dataContext;
      this.selectedCounty = item !== null ? parseInt(item.RS, 10) : -1;
      this.onCountySelected(item);
    });

    this.#map.series.push(newSeries);
    this.#countySeries.set(state.id, newSeries);
    this.#countyHeatLegends.set(state.id, heatLegend);
  }

  /** @private
   *  @param newState {StateMetaData} */
  stateSelected(newState) {
    // If a new state is selected hide the old one.
    if (this.selectedState !== newState.id && this.#countySeries.has(this.selectedState)) {
      this.#countySeries.get(this.selectedState).hide();
      this.#countyHeatLegends.get(this.selectedState).hide();
    }

    if (newState.id > 0) {
      // A state is selected.
      // County series will be loaded lazily.
      if (this.#countySeries.has(newState.id)) {
        this.#countySeries.get(newState.id).show();
        this.#countyHeatLegends.get(newState.id).show();
      } else {
        const state = stateById(newState.id);
        this.createCountySeries(state);
      }

      this.selectedState = newState.id;
      this.#map.zoomToMapObject(this.#stateChoroplethSeries.getPolygonById(newState.id));
      this.#stateHeatLegend.hide();
      this.#stateChoroplethSeries.hide();
      this.#stateBackgroundSeries.show();
    } else {
      // No state is selected.
      this.selectedState = -1;
      this.#map.goHome();
      this.#stateHeatLegend.show();
      this.#stateChoroplethSeries.show();
      this.#stateBackgroundSeries.hide();
    }
  }
}

export class StateMetaData {
  /** @type number */
  id;

  /** @type string | undefined */
  name;

  /** @type string | undefined */
  fileName;
}

/** @type {Array<StateMetaData>} */
const STATE_METADATA = [
  {id: 1, name: 'Schleswig-Holstein', fileName: 'schleswig-holstein'},
  {id: 2, name: 'Hamburg', fileName: 'hamburg'},
  {id: 3, name: 'Niedersachsen', fileName: 'niedersachsen'},
  {id: 4, name: 'Bremen', fileName: 'bremen'},
  {id: 5, name: 'Nordrhein-Westfalen', fileName: 'nordrhein-westfalen'},
  {id: 6, name: 'Hessen', fileName: 'hessen'},
  {id: 7, name: 'Rheinland-Pfalz', fileName: 'rheinland-pfalz'},
  {id: 8, name: 'Baden-Württemberg', fileName: 'baden-wuerttemberg'},
  {id: 9, name: 'Bayern', fileName: 'bayern'},
  {id: 10, name: 'Saarland', fileName: 'saarland'},
  {id: 11, name: 'Berlin', fileName: 'berlin'},
  {id: 12, name: 'Brandenburg', fileName: 'brandenburg'},
  {id: 13, name: 'Mecklenburg-Vorpommern', fileName: 'mecklenburg-vorpommern'},
  {id: 14, name: 'Sachsen', fileName: 'sachsen'},
  {id: 15, name: 'Sachsen-Anhalt', fileName: 'sachsen-anhalt'},
  {id: 16, name: 'Thüringen', fileName: 'thueringen'},
];

/**
 * @param id {number}
 * @return {StateMetaData}
 */
function stateById(id) {
  return STATE_METADATA[id - 1];
}
