import {connect} from "react-redux";

import React from 'react';

import './TimeMap.scss';
import { setSelected } from "../../redux/app";
import InteractiveHeatMap from "../../common/interactive_heat_map";
import {roundToUTCMidnight} from "../../common/utils";

class TimeMap extends React.Component {
  state = {
    /** Map<number, Map<number, number>> */
    stateTimes: new Map(),

    /** Map<number, Map<string, number>> */
    countyTimes: new Map()
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
        if (this.props.selection !== null) {
          this.props.setSelected(null);
        }
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

  /** @private */
  _calcStateData() {
    const times = new Map();

    /** @type Map<number, number> | null */
    let lastStates = null;
    for (let d = this.props.time.startDate; d < this.props.time.endDate; d += 24 * 60 * 60 * 1000) {
      let date = roundToUTCMidnight(d);
      let states = new Map();

      for (let i = 1; i <= 16; i++) {
        const s = this.props.states.all[i];
        const value = s.find(e => e.date >= date);

        if (value) {
          states.set(i, value.Confirmed);
        } else if (lastStates !== null) {
          states.set(i, lastStates.get(i));
        } else {
          states.set(i, 0);
        }
      }

      times.set(date, states);
      lastStates = states;
    }

    if (times.size > 0) {
      this.setState({
        stateTimes: times
      });
    }
  }

  /** @private */
  _calcCountyData() {
    const times = new Map();

    /** @type Map<number, number> | null */
    let lastCounties = null;
    for (let d = this.props.time.startDate; d < this.props.time.endDate; d += 24 * 60 * 60 * 1000) {
      let date = roundToUTCMidnight(d);
      let counties = new Map();

      for (let [id, c] of Object.entries(this.props.counties.all)) {
        const value = c.find(e => e.date >= date);

        if (value) {
          counties.set(id, value.Confirmed);
        } else if (lastCounties !== null) {
          counties.set(id, lastCounties.get(id));
        } else {
          counties.set(id, 0);
        }
      }

      times.set(date, counties);
      lastCounties = counties;
    }

    if (times.size > 0) {
      this.setState({
        countyTimes: times
      });
    }
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.states !== prevProps.states || this.state.stateTimes.size === 0) {
      this._calcStateData();
    }

    if (this.props.counties !== prevProps.counties || this.state.countyTimes.size === 0) {
      this._calcCountyData();
    }

    const currDate = this.props.time.currentDate;
    if (prevProps.time.currentDate !== currDate) {
      if (this.map.selectedState !== -1 && this.state.countyTimes.has(currDate)) {
        this.map.setCountyValues(this.state.countyTimes.get(currDate));
      }
      if (this.state.stateTimes.has(currDate)) {
        this.map.setStateValues(this.state.stateTimes.get(currDate));
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
    time: state.time,
    states: state.app.states,
    counties: state.app.counties,
    populations: state.app.populations.states
  }
}

const ConnectedTimeMap = connect(mapState, { setSelected })(TimeMap);

export {ConnectedTimeMap as TimeMap};