import {createSlice} from '@reduxjs/toolkit';
import {groupBy, renameKey, sumByKey, filterJSObject, stateIdFromCountyId, lastElement} from '../common/utils';
import {setTimeBounds} from './time';

import axios from 'axios';

export const Datasets = {
  STATES: 'states',
  COUNTIES: 'counties',
  GERMANY: 'germany',
};

const dataset2datakey = {
  [Datasets.STATES]: 'ID_State',
  [Datasets.COUNTIES]: 'ID_County',
};

const populationKeyMap = {
  [Datasets.STATES]: 'Bundesland',
  [Datasets.COUNTIES]: 'Landkreis',
};

const slice = createSlice({
  name: 'app',
  initialState: {
    selected: null,
  },
  reducers: {
    init: (state, action) => {
      //return state;
    },
    setSelected: (state, action) => {
      state.selected = action.payload;
    },
  },
});

export const {init, setSelected} = slice.actions;

export const initializeApp = () => (dispatch) => {
  axios
    .get('assets/config.json')
    .then((response) => response.data)
    .then((config) => dispatch(init(config)));
};

export default slice.reducer;
