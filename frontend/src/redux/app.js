import { createSlice } from '@reduxjs/toolkit';
import { groupBy, renameKey, sumByKey } from '../common/utils';
import { setTimeBounds } from "./time";

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
      age: {}
    },
    states: {
      all: {},
      gender: {},
      age: {}
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
      state.selected = action.payload;
      if (action.payload !== null) {
        const population = state.populations[action.payload.dataset].find(
          (e) =>
            e[populationKeyMap[action.payload.dataset]] === action.payload.label
        );

        if (population) {
          state.selected.population = population.EWZ;
        }
      }
    },
    setStateData: (state, action) => {
      const [all, age, gender] = action.payload;
      Object.keys(all).forEach((k) => {
        all[k].sort((a, b) => {
          return a.date - b.date;
        });
      });

      state.states = { all, gender, age };
    },
    setCountyData: (state, action) => {
      const [all, age, gender] = action.payload;
      state.counties = { all, gender, age };
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
};

export const getSelectedData = (state) => {
  const { selected, ...s } = state.app;

  if (selected === null) {
    return null;
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

  const result = {
    all: dataset.all[`${selected.id}`],
    age: dataset.age[`${selected.id}`],
    gender: dataset.gender[`${selected.id}`],
    start: dataset.all[`${selected.id}`][0]
  };

  return result;
};

export default slice.reducer;
