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

  /** @type IHeatRule */
  stateHeatRule = null;

  componentDidMount() {
    this.map = am4core.create("timeMapDiv", am4maps.MapChart);
    this.map.projection = new am4maps.projections.Mercator();

    this.stateSeries = new am4maps.MapPolygonSeries();

    this.stateSeries.geodataSource.url = "assets/german-states.geojson";
    this.stateSeries.useGeodata = true;
    this.map.series.push(this.stateSeries);

    this.stateSeries.dataFields.value = "value";

    this.stateHeatRule = {
      property: "fill",
      target: this.stateSeries.mapPolygons.template,
      min: am4core.color("#EEE", 1),
      max: am4core.color("#F00", 1),
      logarithmic: true,
      dataField: "value"
    }

    this.stateSeries.heatRules.push(this.stateHeatRule);

    const statePolygonTemplate = this.stateSeries.mapPolygons.template;
    statePolygonTemplate.tooltipText = "{name}: {value}";
    // statePolygonTemplate.fill = this.map.colors.getIndex(0);
    statePolygonTemplate.nonScalingStroke = true;

    // Hover state
    // const stateHover = statePolygonTemplate.states.create("hover");
    // stateHover.properties.stroke = am4core.color("#000");

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

    this.map.events.on("zoomlevelchanged", e => {
      if (this.map.zoomLevel === 1) {
        this.props.setSelected(null);
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

        if (this.stateSeries.heatRules.length > 0) {
          this.stateSeries.heatRules.clear();
          this.stateSeries.mapPolygons.template.fill = am4core.color("#FF00FF"); // this.map.colors.getIndex(0);
          // this.chart.deepInvalidate();
          this.stateSeries.deepInvalidate();
        }
      } else {
        this.currentSelected = -1;
        this.map.goHome();

        if (this.stateSeries.heatRules.length === 0) {
          this.stateSeries.heatRules.push(this.stateHeatRule);
        }
      }
  }

  calcStateData() {
    this.state.times = [];

    let lastEntry = null;
    for (let d = this.props.time.startDate; d < this.props.time.endDate; d += 24 * 60 * 60 * 1000) {
      let entry = { day: d, states: {} };

      for (let i = 1; i <= 16; i++) {
        const s = this.props.states.all[i];
        const value = s.find(e => e.date >= d);
        if (value) {
          entry.states[i] = value.Confirmed;
        } else if (lastEntry !== null) {
          entry.states[i] = lastEntry.states[i];
        } else {
          entry.states[i] = 0;
        }
      }

      this.state.times.push(entry);
      lastEntry = entry;
    }
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if ((this.props.states !== prevProps.states || this.state.times.length === 0)) {
      this.calcStateData();
    }

    if (prevProps.time.currentDate !== this.props.time.currentDate) {
      const current = this.state.times.find(t => t.day >= this.props.time.currentDate);
      if (current) {
        for (let s of this.stateSeries.data) {
          s.value = current.states[s.id];
        }
      }

      this.stateSeries.invalidateRawData();
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
    time: state.time,
    states: state.app.states,
    populations: state.app.populations.states
  }
}

const ConnectedTimeMap = connect(mapState, {setSelected})(TimeMap);

export {ConnectedTimeMap as TimeMap};