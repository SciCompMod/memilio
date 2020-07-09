import * as am4core from "@amcharts/amcharts4/core";
import * as am4maps from "@amcharts/amcharts4/maps";


export class InteractiveHeatMap {
  /** @private
   *  @type MapChart */
  _map = null;

  /** @private
   *  @type MapPolygonSeries */
  _stateSeries = null;

  /** @private
   *  @type Map<number, MapPolygonSeries> */
  _countySeries = new Map();

  /** @private
   *  @type IHeatRule */
  _stateHeatRule = null;

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
    this._map = am4core.create(id, am4maps.MapChart);
    this._map.projection = new am4maps.projections.Mercator();
    this._map.hiddenState.transitionDuration = 500;

    this._stateSeries = new am4maps.MapPolygonSeries();

    this._stateSeries.geodataSource.url = "assets/german-states.geojson";
    this._stateSeries.useGeodata = true;
    this._map.series.push(this._stateSeries);

    this._stateHeatRule = {
      property: "fill",
      target: this._stateSeries.mapPolygons.template,
      min: am4core.color("#DDD", 1),
      max: am4core.color("#F00", 1),
      logarithmic: true,
    }

    this._stateSeries.heatRules.push(this._stateHeatRule);

    const statePolygonTemplate = this._stateSeries.mapPolygons.template;
    statePolygonTemplate.tooltipText = "{name}: {value}";
    statePolygonTemplate.nonScalingStroke = true;
    statePolygonTemplate.applyOnClones = true;

    statePolygonTemplate.events.on("hit", e => {
      const item = e.target.dataItem.dataContext;
      this.onStateSelected(item);
      this.onCountySelected(null);
      this.selectedCounty = -1;
      this._stateSelected(item);
    });

    this._map.events.on("zoomlevelchanged", e => {
      if (this._map.zoomLevel === 1) {
        this.onStateSelected(null);
        this.onCountySelected(null);
        this.selectedCounty = -1;
        this._stateSelected({id: -1});
      }
    });
  }

  /** @param values {Map<number, number>} */
  setStateValues(values) {
    for (let stateDatum of this._stateSeries.data) {
      stateDatum.value = values.get(stateDatum.id);
    }

    this._stateSeries.invalidateRawData();
  }

  /** @param values {Map<string, number>} */
  setCountyValues(values) {
    if (this._countySeries.has(this.selectedState)) {
      const series = this._countySeries.get(this.selectedState);
      for (let countyDatum of series.data) {
        countyDatum.value = values.get(parseInt(countyDatum.RS, 10).toString());
      }

      series.invalidateRawData();
    }
  }

  /** @private
   *  @param state {StateMetaData} */
  _createCountySeries(state) {
    const newSeries = new am4maps.MapPolygonSeries();
    newSeries.geodataSource.url = "assets/counties/" + state.fileName + ".geojson";
    newSeries.useGeodata = true;

    const countyPolygonTemplate = newSeries.mapPolygons.template;
    countyPolygonTemplate.tooltipText = "{id}: {value}";
    countyPolygonTemplate.nonScalingStroke = true;
    countyPolygonTemplate.applyOnClones = true;

    const countyHeatRule = {
      property: "fill",
      target: countyPolygonTemplate,
      min: am4core.color("#DDD", 1),
      max: am4core.color("#F00", 1),
      logarithmic: true,
    }

    newSeries.heatRules.push(countyHeatRule);

    countyPolygonTemplate.events.on("hit", e => {
      const item = e.target.dataItem.dataContext;
      this.selectedCounty = item !== null ? parseInt(item.RS, 10) : -1;
      this.onCountySelected(item);
    });

    this._map.series.push(newSeries);
    this._countySeries.set(state.id, newSeries);
  }

  /** @private
   *  @param newState {StateMetaData} */
  _stateSelected(newState) {
    // If a new state is selected hide the old one.
    if (this.selectedState !== newState.id && this._countySeries.has(this.selectedState)) {
      this._countySeries.get(this.selectedState).hide();
    }

    if (newState.id > 0) { // A state is selected.
      // County series will be loaded lazily.
      if (this._countySeries.has(newState.id)) {
        this._countySeries.get(newState.id).show();
      } else {
        const state = stateById(newState.id);
        this._createCountySeries(state);
      }

      this.selectedState = newState.id;
      this._map.zoomToMapObject(this._stateSeries.getPolygonById(newState.id));
      this._stateSeries.hide();
    } else { // No state is selected.
      this.selectedState = -1;
      this._map.goHome();
      this._stateSeries.show();
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
  {id: 1, name: "Schleswig-Holstein", fileName: "schleswig-holstein"},
  {id: 2, name: "Hamburg", fileName: "hamburg"},
  {id: 3, name: "Niedersachsen", fileName: "niedersachsen"},
  {id: 4, name: "Bremen", fileName: "bremen"},
  {id: 5, name: "Nordrhein-Westfalen", fileName: "nordrhein-westfalen"},
  {id: 6, name: "Hessen", fileName: "hessen"},
  {id: 7, name: "Rheinland-Pfalz", fileName: "rheinland-pfalz"},
  {id: 8, name: "Baden-Württemberg", fileName: "baden-wuerttemberg"},
  {id: 9, name: "Bayern", fileName: "bayern"},
  {id: 10, name: "Saarland", fileName: "saarland"},
  {id: 11, name: "Berlin", fileName: "berlin"},
  {id: 12, name: "Brandenburg", fileName: "brandenburg"},
  {id: 13, name: "Mecklenburg-Vorpommern", fileName: "mecklenburg-vorpommern"},
  {id: 14, name: "Sachsen", fileName: "sachsen"},
  {id: 15, name: "Sachsen-Anhalt", fileName: "sachsen-anhalt"},
  {id: 16, name: "Thüringen", fileName: "thueringen"}
]

/**
 * @param id {number}
 * @return {StateMetaData}
 */
function stateById(id) {
  return STATE_METADATA[id - 1];
}