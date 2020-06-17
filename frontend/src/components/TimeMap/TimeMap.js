import {connect} from "react-redux";

import React from 'react';
import * as am4core from "@amcharts/amcharts4/core";
import * as am4maps from "@amcharts/amcharts4/maps";

import './TimeMap.scss';
import {setSelected} from "../../redux/app";

class StateMetaData {
  id; name; fileName;
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

class TimeMap extends React.Component {
  state = {};

  /** @type MapChart */
  map = null;

  /** @type MapPolygonSeries */
  stateSeries = null;

  /** @type Map<number, MapPolygonSeries> */
  countySeries = new Map();

  /** @type number */
  currentSelected = -1;

  componentDidMount() {
    this.map = am4core.create("timeMapDiv", am4maps.MapChart);
    this.map.projection = new am4maps.projections.Mercator();

    this.stateSeries = new am4maps.MapPolygonSeries();
    this.stateSeries.geodataSource.url = "assets/german-states.geojson";
    this.stateSeries.useGeodata = true;
    this.map.series.push(this.stateSeries);

    const statePolygonTemplate = this.stateSeries.mapPolygons.template;
    statePolygonTemplate.tooltipText = "{name}";
    statePolygonTemplate.fill = this.map.colors.getIndex(0);
    statePolygonTemplate.nonScalingStroke = true;

    // Hover state
    const stateHover = statePolygonTemplate.states.create("hover");
    stateHover.properties.fill = am4core.color("#367B25");

    statePolygonTemplate.events.on("hit", ev => {
      const item = ev.target.dataItem.dataContext;
      this.stateSelected(item);
      this.props.setSelected({
        dataset: "states",
        id: item.id,
        label: item.name,
        population: item.destatis.population
      });
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

    countyPolygonTemplate.events.on("hit", ev => {
      const item = ev.target.dataItem.dataContext;
      this.props.setSelected({
        dataset: "counties",
        id: parseInt(item.RS, 10),
        label: item.id,
        population: item.destatis.population
      });
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
      } else {
        this.currentSelected = -1;
        this.map.goHome();
      }
  }


  render(ctx) {
    return (
      <div style={{height: "100%"}}>
        <div id="timeMapDiv" style={{height: "100%"}}/>
      </div>
    )
  }
}

function mapState(state) {
  return {
    selection: state.app.selected,
    time: state.time
  }
}

const ConnectedTimeMap = connect(mapState, {setSelected})(TimeMap);

export {ConnectedTimeMap as TimeMap};