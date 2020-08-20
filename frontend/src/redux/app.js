import {createSlice} from '@reduxjs/toolkit';
import axios from 'axios';

export const Datasets = {
  STATES: 'states',
  COUNTIES: 'counties',
  GERMANY: 'germany',
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
