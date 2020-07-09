import {connect} from "react-redux";

import React from 'react';

import './TimeMap.scss';
import {setSelected} from "../../redux/app";
import {InteractiveHeatMap} from "../../common/interactive_heat_map.js";

class TimeMap extends React.Component {
  state = {
    times: []
  };

  /** @type InteractiveHeatMap */
  map = null;

  componentDidMount() {
    this.map = new InteractiveHeatMap("timeMapDiv");

    this.map.onStateSelected = newState => {
      if (newState !== null) {
        this.props.setSelected({
          dataset: "states",
          id: newState.id,
          label: newState.name,
          population: newState.destatis.population
        });
      } else {
        this.props.setSelected(null);
      }
    };

    this.map.onCountySelected = newCounty => {
      if (newCounty !== null) {
        this.props.setSelected({
          dataset: "counties",
          id: parseInt(newCounty.RS, 10),
          label: newCounty.id,
          population: newCounty.destatis.population
        });
      }
    };
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
        for (let s of this.map.stateSeries.data) {
          s.value = current.states[s.id];
        }
      }

      this.map.stateSeries.invalidateRawData();
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