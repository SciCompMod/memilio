import * as am4core from "@amcharts/amcharts4/core";
import * as am4maps from "@amcharts/amcharts4/maps";


export class InteractiveHeatMap {
  /** @type MapChart */
  map = null;

  /** @type MapPolygonSeries */
  stateSeries = null;

  /** @type Map<number, MapPolygonSeries> */
  countySeries = new Map();

  /** @type number */
  currentSelected = -1;

  /** @type IHeatRule */
  stateHeatRule = null;

  /**
   * @callback StateSelectedCallback
   * @param newState {StateMetaData}
   */
  onStateSelected = (newState) => {};

  /**
   * @callback CountySelectedCallback
   * @param newCounty {any}
   */
  onCountySelected = (newCounty) => {};

  constructor(id) {
    this.map = am4core.create(id, am4maps.MapChart);
    this.map.projection = new am4maps.projections.Mercator();
    this.map.hiddenState.transitionDuration = 500;

    this.stateSeries = new am4maps.MapPolygonSeries();

    this.stateSeries.geodataSource.url = "assets/german-states.geojson";
    this.stateSeries.useGeodata = true;
    this.map.series.push(this.stateSeries);

    this.stateHeatRule = {
      property: "fill",
      target: this.stateSeries.mapPolygons.template,
      min: am4core.color("#DDD", 1),
      max: am4core.color("#F00", 1),
      logarithmic: true,
    }

    this.stateSeries.heatRules.push(this.stateHeatRule);

    const statePolygonTemplate = this.stateSeries.mapPolygons.template;
    statePolygonTemplate.tooltipText = "{name}: {value} Fälle";
    statePolygonTemplate.nonScalingStroke = true;
    statePolygonTemplate.applyOnClones = true;

    statePolygonTemplate.events.on("hit", e => {
      const item = e.target.dataItem.dataContext;
      this.stateSelected(item);
      this.onStateSelected(item);
      this.onCountySelected(null);
    });

    this.map.events.on("zoomlevelchanged", e => {
      if (this.map.zoomLevel === 1) {
        this.onStateSelected(null);
        this.onCountySelected(null);
        this.stateSelected({id: -1});
      }
    });
  }

  /** @param state {StateMetaData} */
  createCountySeries(state) {
    const newSeries = new am4maps.MapPolygonSeries();
    newSeries.geodataSource.url = "assets/counties/" + state.fileName + ".geojson";
    newSeries.useGeodata = true;

    const countyPolygonTemplate = newSeries.mapPolygons.template;
    countyPolygonTemplate.tooltipText = "{id}";
    countyPolygonTemplate.fill = this.map.colors.getIndex(0);
    countyPolygonTemplate.nonScalingStroke = true;

    const hs = countyPolygonTemplate.states.create("hover");
    hs.properties.fill = am4core.color("#367B25");

    countyPolygonTemplate.events.on("hit", e => {
      const item = e.target.dataItem.dataContext;
      this.onCountySelected(item);
    });

    this.map.series.push(newSeries);
    this.countySeries.set(state.id, newSeries);
  }

  /** @param newState {StateMetaData} */
  stateSelected(newState) {
    if (this.currentSelected !== newState.id) {
      if (this.countySeries.has(this.currentSelected)) {
        this.countySeries.get(this.currentSelected).hide();
      }
    }

    if (newState.id > 0) {
      if (this.countySeries.has(newState.id)) {
        this.countySeries.get(newState.id).show();
      } else {
        const state = stateById(newState.id);
        this.createCountySeries(state);
      }
      this.currentSelected = newState.id;
      this.map.zoomToMapObject(this.stateSeries.getPolygonById(newState.id));

      if (this.stateSeries.heatRules.length > 0) {
        this.stateSeries.heatRules.clear();
        this.stateSeries.mapPolygons.template.fill = am4core.color("#EEE");
      }
    } else {
      this.currentSelected = -1;
      this.map.goHome();

      if (this.stateSeries.heatRules.length === 0) {
        this.stateSeries.heatRules.push(this.stateHeatRule);
      }
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