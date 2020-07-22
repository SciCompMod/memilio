import {createSlice} from '@reduxjs/toolkit';
import {groupBy, renameKey, sumByKey, filterJSObject, stateIdFromCountyId} from '../common/utils';
import {setTimeBounds} from "./time";

import axios from 'axios';

export const Datasets = {
  STATES: 'states',
  COUNTIES: 'counties',
  GERMANY: 'germany'
};

const dataset2datakey = {
  [Datasets.STATES]: 'ID_State',
  [Datasets.COUNTIES]: 'ID_County'
};

const populationKeyMap = {
  [Datasets.STATES]: 'Bundesland',
  [Datasets.COUNTIES]: 'Landkreis'
};

const slice = createSlice({
  name: 'app',
  initialState: {
    selected: null,
    populations: {
      states: [],
      counties: [],
      total: 0
    },
    germany: {
      all: {},
      gender: {},
      age: {},
    },
    states: {
      all: {},
      gender: {},
      age: {},
    },
    counties: {
      all: {},
      gender: {},
      age: {}
    }
  },
  reducers: {
    init: (state, action) => {
      //return state;
    },
    setSelected: (state, action) => {
      if (action.payload !== null) {
        const population = state.populations[action.payload.dataset].find(
          (e) =>
            e[populationKeyMap[action.payload.dataset]] === action.payload.label
        );

        if (population) {
          state.selected.population = population.EWZ;
        }
      } else {
        state.selected = {dataset: "germany", id: 0, label: "Germany", population: state.populations.total};
      }
    },
    setCountryData(state, action) {
      state.germany = action.payload;
    },
    setStateData: (state, action) => {
      const [all, age, gender] = action.payload;
      Object.keys(all).forEach((k) => {
        all[k].sort((a, b) => {
          return a.date - b.date;
        });
      });

      state.states = {all, gender, age};
    },
    setCountyData: (state, action) => {
      const [all, age, gender] = action.payload;
      state.counties = {all, gender, age};
    },
    setPopulations: (state, action) => {
      state.populations = action.payload;
      state.populations.total = sumByKey(state.populations.states, 'EWZ');
    }
  }
});

export const {
  init,
  setSelected,
  setCountryData,
  setStateData,
  setCountyData,
  setPopulations
} = slice.actions;

export const initializeApp = () => (dispatch) => {
  axios
    .get('assets/config.json')
    .then((response) => response.data)
    .then((config) => dispatch(init(config)));
};

export const fetchData = () => async (dispatch) => {
  // load state data
  let state = await Promise.all([
    axios.get('assets/all_state_rki.json').then((response) => response.data),
    axios
      .get('assets/all_state_age_rki.json')
      .then((response) => response.data),
    axios
      .get('assets/all_state_gender_rki.json')
      .then((response) => response.data)
  ]);

  state = state.map((s) => {
    return groupBy(
      renameKey(s, 'Date', 'date'),
      dataset2datakey[Datasets.STATES]
    );
  });

  dispatch(setStateData(state));

  // Generate Germany data from state data by summing it up.
  function sumByState(stateData) {
    const summedData = [];
    for (let stateEntry of Object.values(stateData)) {
      for (let timeStamp of stateEntry) {
        const result = summedData.find(entry => entry.date === timeStamp.date);
        if (result) {
          result.Confirmed += timeStamp.Confirmed;
          result.Deaths += timeStamp.Deaths;
          result.Recovered += timeStamp.Recovered;
        } else {
          summedData.push({
            ID_Country: 0,
            Country: "Germany",
            Confirmed: timeStamp.Confirmed,
            Deaths: timeStamp.Deaths,
            Recovered: timeStamp.Recovered,
            date: timeStamp.date
          });
        }
      }
    }
    summedData.sort((a, b) => a.date - b.date);
    return summedData;
  }

  const [all, age, gender] = state;

  dispatch(setCountryData({all: sumByState(all), gender: sumByState(gender), age: sumByState(age)}));

  // TODO This is a temporary hack!
  dispatch(setTimeBounds({start: state[0][9][0].date, end: state[0][9][119].date}))

  // load county data
  let county = await Promise.all([
    axios.get('assets/all_county_rki.json').then((response) => response.data),
    axios
      .get('assets/all_county_age_rki.json')
      .then((response) => response.data),
    axios
      .get('assets/all_county_gender_rki.json')
      .then((response) => response.data)
  ]);

  county = county.map((s) => {
    return groupBy(
      renameKey(s, 'Date', 'date'),
      dataset2datakey[Datasets.COUNTIES]
    );
  });

  dispatch(setCountyData(county));

  const populations = await axios
    .get('assets/populations.json')
    .then((response) => response.data);

  dispatch(setPopulations(populations));
  dispatch(setSelected(null));
};

export const getSelectedData = (state) => {
  const {selected, ...s} = state.app;

  if (selected === null) {
    return null;
  }

  if (selected.dataset === "germany") {
    return {...s.germany, start: s.germany.all[0]};
  }

  const dataset = s[selected.dataset];

  if (!dataset) {
    return null;
  }

  if (
    ![dataset.all, dataset.age, dataset.gender].every((d) =>
      d.hasOwnProperty(`${selected.id}`)
    )
  ) {
    return null;
  }

  return {
    all: dataset.all[`${selected.id}`],
    age: dataset.age[`${selected.id}`],
    gender: dataset.gender[`${selected.id}`],
    start: dataset.all[`${selected.id}`][0]
  };
};

/**
 * Returns the RKI data of the selected regions children.
 *
 * @param state The application state.
 * @return {{all: any age: any, gender: any}|null}
 */
export const getSelectedChildData = state => {
  const {selected, ...s} = state.app;

  if (selected === null) {
    return null;
  }

  if (selected.dataset === "germany") {
    const dataset = s.states;
    return {
      all: dataset.all,
      age: dataset.age,
      gender: dataset.gender
    };
  }

  /**
   * County ids are defined as state id (SS) and three numbers (CCC) => SSCCC. We can check if a county belongs to a
   * state by comparing the state id to the first two numbers of the county id.
   * @param stateId {number}
   * @param counties {{all: any age: any, gender: any}}
   * @return {{all: any age: any, gender: any}}
   */
  const filterByStateId = (stateId, counties) => filterJSObject(counties, (id, _) => stateIdFromCountyId(id) === stateId);

  if (selected.dataset === "states" || selected.dataset === "counties") {
    const stateId = selected.dataset === "states" ? selected.id : stateIdFromCountyId(selected.id);
    const dataset = s.counties;
    return {
      all: filterByStateId(stateId, dataset.all),
      age: filterByStateId(stateId, dataset.age),
      gender: filterByStateId(stateId, dataset.gender),
    };
  }

  return null;
};

/**
 * Returns the populations of a regions child regions. E.g. If the regionId is 0 (Germany) the populations of all states
 * are returned. Possible cases are:
 *
 * regionId belongs to germany  => state ids
 * regionId belongs to a state  => county ids of that state
 * regionId belongs to a county => county ids of counties in the same state
 *
 * @param state Application state.
 * @param regionId {number} A state or county id.
 * @return Map<number, number>
 */
export function getPopulationsOfRegion(state, regionId) {
  if (regionId === 0) {
    const result = new Map();
    for (let stateRegion of state.app.populations.states) {
      result.set(stateRegion.stateKey, stateRegion.EWZ);
    }
    return result;
  }

  let counties;
  if (regionId < 100) {
    counties = state.app.populations.counties.filter(county => stateIdFromCountyId(county.countyKey) === regionId);
  } else {
    counties = state.app.populations.counties.filter(county => stateIdFromCountyId(county.countyKey) === stateIdFromCountyId(regionId));
  }

  const result = new Map();
  for (let county of counties) {
    result.set(county.countyKey, county.EWZ);
  }

  return result;
}

export default slice.reducer;
