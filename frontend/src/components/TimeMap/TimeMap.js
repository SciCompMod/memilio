import {connect} from 'react-redux';

import React from 'react';

import {setSelected} from '../../redux/app';
import InteractiveHeatMap from '../../common/interactive_heat_map';
import {roundToUTCMidnight} from '../../common/utils';
import {RKIDatastore as rki, Tables} from '../../common/rki-datastore';

import './TimeMap.scss';
/**
 * This Component has two major functions:
 * 1. Select different Regions for the chart.
 * 2. Display a choropleth map of the infection data over time.
 *    Currently the choropleth map displays either the RKI data, or if available simulated SEIR data of all infected.
 *    A selection mechanism will be integrated at a later time.
 *    If the timeline changes the current time the map will update with the values of the selected date.
 */
class TimeMap extends React.Component {
  state = {
    /** @type Map<number, Map<number, number>> */
    stateTimes: new Map(),

    /** @type Map<number, Map<string, number>> */
    countyTimes: new Map(),

    /** @type Map<number, Map<string, number>> */
    seirTimes: new Map(),
  };

  /** @type InteractiveHeatMap */
  map = null;

  componentDidMount() {
    this.map = new InteractiveHeatMap('timeMapDiv');

    this.map.onStateSelected = (newState) => {
      if (newState !== null) {
        this.props.setSelected({
          dataset: 'states',
          id: newState.id,
          label: newState.name,
          population: newState.destatis.population,
        });
      } else {
        if (this.props.selection !== null) {
          //this.props.setSelected(null);
        }
      }
    };

    this.map.onCountySelected = (newCounty) => {
      if (newCounty !== null) {
        this.props.setSelected({
          dataset: 'counties',
          id: parseInt(newCounty.RS, 10),
          label: newCounty.id,
          population: newCounty.destatis.population,
        });
      }
    };
  }

  /** @private */
  async _calcStateData() {
    const times = new Map();

    /** @type Map<number, number> | null */
    let lastStates = null;
    for (let d = this.props.time.startDate; d < this.props.time.endDate; d += 24 * 60 * 60 * 1000) {
      let date = roundToUTCMidnight(d);
      let states = new Map();

      const s = await rki.get(Tables.STATES, {start: date, end: date});
      for (let stateID = 1; stateID <= 16; stateID++) {
        const value = s.find((e) => e.ID_State === stateID);
        if (value) {
          states.set(stateID, value.Confirmed - value.Recovered);
        } else if (lastStates !== null) {
          states.set(stateID, lastStates.get(stateID));
        } // else {
        //  states.set(stateID, 0);
        //}
      }

      //for (let i = 1; i <= 16; i++) {
      // const s = this.props.states.all[i];
      //const value = s.find((e) => e.date >= date);
      //const s = await rki.getState(i, {start: date, end: date});
      //console.log(s);
      //continue;
      // const value =
      /*if (value) {
          states.set(i, value.Confirmed);
        } else if (lastStates !== null) {
          states.set(i, lastStates.get(i));
        } else {
          states.set(i, 0);
        }*/
      //}

      times.set(date, states);
      lastStates = states;
    }

    if (times.size > 0) {
      this.setState({
        stateTimes: times,
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
        const value = c.find((e) => e.date >= date);

        if (value) {
          counties.set(parseInt(id), value.Confirmed);
        } else if (lastCounties !== null) {
          counties.set(parseInt(id), lastCounties.get(parseInt(id)));
        } else {
          counties.set(parseInt(id), 0);
        }
      }

      times.set(date, counties);
      lastCounties = counties;
    }

    if (times.size > 0) {
      this.setState({
        countyTimes: times,
      });
    }
  }

  /** @private */
  _calcSeirData() {
    const times = new Map();

    /** @type Map<number, number> | null */
    let lastRegions = null;
    for (let d = this.props.time.startDate; d < this.props.time.endDate; d += 24 * 60 * 60 * 1000) {
      let date = roundToUTCMidnight(d);
      let regions = new Map();

      for (let [id, region] of Object.entries(this.props.seirRegions)) {
        const value = region.find((e) => e.date >= date);

        if (value) {
          regions.set(parseInt(id), value.I);
        } else if (lastRegions !== null) {
          regions.set(parseInt(id), lastRegions.get(parseInt(id)));
        } else {
          regions.set(parseInt(id), 0);
        }
      }

      times.set(date, regions);
      lastRegions = regions;
    }

    if (times.size > 0) {
      this.setState({
        seirTimes: times,
      });
    }
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (
      this.props.states !== prevProps.states ||
      this.state.stateTimes.size === 0 ||
      this.props.time.startDate !== prevProps.time.startDate ||
      this.props.time.endDate !== prevProps.time.endDate
    ) {
      this._calcStateData();
    }

    if (
      this.props.counties !== prevProps.counties ||
      this.state.countyTimes.size === 0 ||
      this.props.time.startDate !== prevProps.time.startDate ||
      this.props.time.endDate !== prevProps.time.endDate
    ) {
      //this._calcCountyData();
    }

    if (
      this.props.seirRegions !== null &&
      (this.props.seirRegions !== prevProps.seirRegions || this.state.seirTimes.size === 0)
    ) {
      //this._calcSeirData();
    }

    const currDate = this.props.time.currentDate;
    if (prevProps.time.currentDate !== currDate) {
      if (this.props.seirRegions !== null) {
        this.map.setDataSetName('SEIR');
        if (this.map.selectedState !== -1) {
          this.map.setCountyValues(this.state.seirTimes.get(currDate));
        } else {
          this.map.setStateValues(this.state.seirTimes.get(currDate));
        }
      } else {
        this.map.setDataSetName('RKI');
        if (this.map.selectedState !== -1 && this.state.countyTimes.has(currDate)) {
          this.map.setCountyValues(this.state.countyTimes.get(currDate));
        }

        if (this.state.stateTimes.has(currDate)) {
          this.map.setStateValues(this.state.stateTimes.get(currDate));
        }
      }
    }
  }

  render(ctx) {
    return (
      <div style={{height: '100%'}}>
        <div id="timeMapDiv" style={{height: '100%'}} />
      </div>
    );
  }
}

function mapState(state) {
  return {
    selection: state.app.selected,
    time: state.time,
    seirRegions: state.seir.regionData,
  };
}

const ConnectedTimeMap = connect(mapState, {setSelected})(TimeMap);

export {ConnectedTimeMap as TimeMap};
