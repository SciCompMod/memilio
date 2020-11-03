import {connect} from 'react-redux';

import React from 'react';

import {setSelected} from '../../redux/app';
import InteractiveHeatMap from '../../common/interactive_heat_map';
import {roundToUTCMidnight} from '../../common/utils';
import rki from '../../common/datastore/sql/rki-sql-store';

import './TimeMap.scss';
import {ButtonGroup} from 'reactstrap';
import {Button} from 'reactstrap';
import {withTranslation} from 'react-i18next';

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
          label: `${newCounty.BEZ} ${newCounty.GEN}`,
          population: newCounty.destatis.population,
        });
      }
    };

    this.setState({dataset: 'states'});
  }

  /** @private */
  async calcStateData(data) {
    const times = new Map();

    let currDate = -1;
    let states = null;

    for (const {stateId, confirmed, recovered, date} of data) {
      const d = roundToUTCMidnight(date);

      if (d !== currDate) {
        times.set(currDate, states);
        currDate = roundToUTCMidnight(date);
        states = new Map(states); // clone previous date so there are no holes in data
      }

      if (stateId > 0) {
        states.set(stateId, confirmed - recovered);
      }
    }

    if (times.size > 0) {
      this.setState({
        stateTimes: times,
      });
    }
  }

  /** @private */
  async calcCountyData(data) {
    const times = new Map();

    let currDate = -1;
    let counties = null;
    for (const {countyId, confirmed, recovered, date} of data) {
      const d = roundToUTCMidnight(date);
      if (d !== currDate) {
        times.set(currDate, counties);
        currDate = roundToUTCMidnight(date);
        counties = new Map(counties); // clone previous date so there are no holes in data
      }

      counties.set(countyId, confirmed - recovered);
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
      this.state.dataset === 'states' &&
      (this.state.dataset !== prevState.dataset ||
        this.props.data.states !== prevProps.data.states ||
        this.state.stateTimes.size === 0 ||
        this.props.time.startDate !== prevProps.time.startDate ||
        this.props.time.endDate !== prevProps.time.endDate)
    ) {
      this.props.data.states.then((data) => this.calcStateData(data));
      this.#map.showStateSeries();
    }

    if (
      this.state.dataset === 'counties' &&
      (this.state.dataset !== prevState.dataset ||
        this.props.data.counties !== prevProps.data.counties ||
        this.state.countyTimes.size === 0 ||
        this.props.time.startDate !== prevProps.time.startDate ||
        this.props.time.endDate !== prevProps.time.endDate)
    ) {
      this.props.data.counties.then((data) => this.calcCountyData(data));
      this.#map.showCountySeries();
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
        if (this.#map.selectedState === -1 && this.state.countyTimes.has(currDate)) {
          this.#map.setCountyValues(this.state.countyTimes.get(currDate));
        }

        if (this.state.stateTimes.has(currDate)) {
          this.#map.setStateValues(this.state.stateTimes.get(currDate));
        }
      }
    }
  }

  render(ctx) {
    const {t} = this.props;
    return (
      <div style={{height: '100%'}}>
        <div className="container-fluid h-100" style={{marginTop: '5px'}}>
          <div className="row justify-content-center h-100">
            <div className="col-12">
              <div className="h-100 d-flex flex-column">
                <div className="row">
                  <ButtonGroup className="col-8 offset-2">
                    <Button
                      color="primary"
                      onClick={() => this.setState({dataset: 'states'})}
                      active={this.state.dataset === 'states'}
                      style={{width: '50%'}}
                    >
                      {t('states')}
                    </Button>
                    <Button
                      color="primary"
                      onClick={() => this.setState({dataset: 'counties'})}
                      active={this.state.dataset === 'counties'}
                      style={{width: '50%'}}
                    >
                      {t('counties')}
                    </Button>
                  </ButtonGroup>
                </div>
                <div className="row flex-grow-1">
                  <div id="timeMapDiv" className="col-12" />
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    );
  }
}

function mapState(state) {
  return {
    data: {
      states: rki.getAllStatesInRange({start: state.time.startDate, end: state.time.endDate}),
      counties: rki.getAllCountiesInRange({start: state.time.startDate, end: state.time.endDate}),
    },
    selection: state.app.selected,
    time: state.time,
    seirRegions: state.seir.regionData,
  };
}

const ConnectedTimeMap = connect(mapState, {setSelected})(TimeMap);
const Translated = withTranslation()(ConnectedTimeMap);

export {Translated as TimeMap};
