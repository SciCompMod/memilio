import {connect} from "react-redux";

import React from 'react';
import * as am4core from "@amcharts/amcharts4/core";
import * as am4maps from "@amcharts/amcharts4/maps";
import am4geodata_gerLow from "@amcharts/amcharts4-geodata/germanyLow";

import './TimeMap.scss';

const STATE_METADATA = {

}

class TimeMap extends React.Component {
  state = {

  };

  /** @type MapChart */
  map = null;

  /** @type MapPolygonSeries */
  stateSeries = null;

  /** @type MapPolygonSeries */
  countySeries = null;

  componentDidMount() {
    this.map = am4core.create("timeMapDiv", am4maps.MapChart);
    this.map.projection = new am4maps.projections.Mercator();

    this.stateSeries = new am4maps.MapPolygonSeries();
    this.stateSeries.geodata = am4geodata_gerLow;
    this.stateSeries.useGeodata = true;
    this.map.series.push(this.stateSeries);

    const statePolygonTemplate = this.stateSeries.mapPolygons.template;
    statePolygonTemplate.tooltipText = "{name}";
    statePolygonTemplate.fill = this.map.colors.getIndex(0);
    statePolygonTemplate.nonScalingStroke = true;

    // Hover state
    const stateHover = statePolygonTemplate.states.create("hover");
    stateHover.properties.fill = am4core.color("#367B25");

    this.countySeries = new am4maps.MapPolygonSeries();
    this.countySeries.geodataSource.url = "assets/landkreise.geojson"
    this.countySeries.useGeodata = true;
    this.countySeries.hidden = true;
    this.map.series.push(this.countySeries);

    const countyPolygonTemplate = this.countySeries.mapPolygons.template;
    countyPolygonTemplate.tooltipText = "{id}";
    countyPolygonTemplate.fill = this.map.colors.getIndex(0);
    countyPolygonTemplate.nonScalingStroke = true;

    // Hover state
    const hs = countyPolygonTemplate.states.create("hover");
    hs.properties.fill = am4core.color("#367B25");
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.selection !== prevProps.selection) {
      if (this.props.selection === null) {
        this.countySeries.hide();
      } else {
        this.countySeries.show();
      }
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

const ConnectedTimeMap = connect(mapState, null)(TimeMap);

export {ConnectedTimeMap as TimeMap};