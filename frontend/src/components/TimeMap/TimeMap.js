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
  #map = null;

  componentDidMount() {
    this.#map = new InteractiveHeatMap('timeMapDiv');
    this.#map.onStateSelected = (newState) => {
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

    this.#map.onCountySelected = (newCounty) => {
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
  async calcStateData() {
    const times = new Map();

    const result = await rki.get(Tables.STATES, {
      start: this.props.time.startDate,
      end: this.props.time.endDate,
      sort: 'date',
    });

    let currDate = -1;
    let states = null;
    let lastStates = null;
    for (const {ID_State, Confirmed, Recovered, date} of result) {
      const d = roundToUTCMidnight(date);
      if (d !== currDate) {
        times.set(currDate, states);
        currDate = roundToUTCMidnight(date);

        lastStates = states;
        states = new Map();
      }

      if (ID_State > 0) {
        states.set(ID_State, Confirmed - Recovered);
      }
    }

    if (times.size > 0) {
      this.setState({
        stateTimes: times,
      });
    }
  }

  /** @private */
  async calcCountyData() {
    const times = new Map();

    const result = await rki.get(Tables.COUNTIES, {
      start: this.props.time.startDate,
      end: this.props.time.endDate,
      sort: 'date',
    });

    let currDate = -1;
    let counties = null;
    let lastCounties = null;
    for (const {ID_County, Confirmed, Recovered, date} of result) {
      const d = roundToUTCMidnight(date);
      if (d !== currDate) {
        times.set(currDate, counties);
        currDate = roundToUTCMidnight(date);

        lastCounties = counties;
        counties = new Map();
      }

      counties.set(ID_County, Confirmed - Recovered);
    }

    if (times.size > 0) {
      this.setState({
        countyTimes: times,
      });
    }
  }

  /** @private */
  calcSeirData() {
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
      this.calcStateData();
    }

    if (
      this.props.counties !== prevProps.counties ||
      this.state.countyTimes.size === 0 ||
      this.props.time.startDate !== prevProps.time.startDate ||
      this.props.time.endDate !== prevProps.time.endDate
    ) {
      this.calcCountyData();
    }

    if (
      this.props.seirRegions !== null &&
      (this.props.seirRegions !== prevProps.seirRegions || this.state.seirTimes.size === 0)
    ) {
      //this.calcSeirData();
    }

    const currDate = this.props.time.currentDate;
    if (prevProps.time.currentDate !== currDate) {
      if (this.props.seirRegions !== null) {
        this.#map.setDataSetName('SEIR');
        if (this.#map.selectedState !== -1) {
          this.#map.setCountyValues(this.state.seirTimes.get(currDate));
        } else {
          this.#map.setStateValues(this.state.seirTimes.get(currDate));
        }
      } else {
        this.#map.setDataSetName('RKI');
        if (this.#map.selectedState !== -1 && this.state.countyTimes.has(currDate)) {
          this.#map.setCountyValues(this.state.countyTimes.get(currDate));
        }

        if (this.state.stateTimes.has(currDate)) {
          this.#map.setStateValues(this.state.stateTimes.get(currDate));
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
    counties: null,
    selection: state.app.selected,
    time: state.time,
    seirRegions: state.seir.regionData,
  };
}

const ConnectedTimeMap = connect(mapState, {setSelected})(TimeMap);

export {ConnectedTimeMap as TimeMap};
