import { createSlice } from '@reduxjs/toolkit';
import { groupBy } from '../common/utils';

import axios from 'axios';

const slice = createSlice({
  name: 'app',
  initialState: {
    selected: {
      dataset: 'states.all',
      entry: 'Bayern'
    },
    states: {
      all: [],
      gender: [],
      age: []
    },
    counties: {
      all: [],
      gender: [],
      age: []
    }
  },
  reducers: {
    init: (state, action) => {
      //return state;
    },
    setStateData: (state, action) => {
      const [all, age, gender] = action.payload;
      state.states = { all, gender, age };
      //console.log(groupBy(age, 'Altersgruppe'));
    },
    setCountyData: (state, action) => {
      const [all, age, gender] = action.payload;
      state.counties = { all, gender, age };
    }
  }
});

export const { init, setStateData, setCountyData } = slice.actions;

export const initializeApp = () => (dispatch) => {
  axios
    .get('/assets/config.json')
    .then((response) => response.data)
    .then((config) => dispatch(init(config)));
};

export const fetchData = () => async (dispatch) => {
  // load state data
  const state = await Promise.all([
    axios.get('/assets/all_state.json').then((response) => response.data),
    axios.get('/assets/all_state_age.json').then((response) => response.data),
    axios.get('/assets/all_state_gender.json').then((response) => response.data)
  ]);

  dispatch(setStateData(state));

  // load county data
  const county = await Promise.all([
    axios.get('/assets/all_county.json').then((response) => response.data),
    axios.get('/assets/all_county_age.json').then((response) => response.data),
    axios
      .get('/assets/all_county_gender.json')
      .then((response) => response.data)
  ]);

  dispatch(setCountyData(county));
};

export const getSelectedData = (state) => {
  const { selected, ...s } = state.app;
  const [a, b] = selected.dataset.split('.');

  return s[a][b].filter((e) => selected.entry == e.Bundesland);
};

export default slice.reducer;
