import { createSlice } from '@reduxjs/toolkit';

const slice = createSlice({
  name: 'app',
  initialState: {},
  reducers: {
    init: (state, action) => {
      return state;
    }
  }
});

export const { init } = slice.actions;

export const initializeState = () => (dispatch) => {
  setTimeout(() => {
    dispatch(
      init({
        measures: [
          {
            label: 'measures.homeoffice',
            damping: 0.7,
            intervals: [
              {
                start: 1586988000000,
                end: 1587679200000
              }
            ]
          },
          {
            label: 'measures.schoolclosing',
            damping: 0.8
          },
          {
            label: 'measures.socialdistancing',
            damping: 0.4
          }
        ],
        parameters: [
          {
            label: 'parameters.contactrate',
            name: 'contact_rate',
            type: 'slider',
            min: 0,
            max: 10,
            step: 0.1,
            default: 2.7
          },
          {
            label: 'parameters.incubation',
            name: 'incubation',
            type: 'slider',
            min: 0,
            max: 10,
            step: 0.1,
            default: 5.2
          },
          {
            label: 'parameters.infection',
            name: 'infection',
            type: 'slider',
            min: 0,
            max: 10,
            step: 0.1,
            default: 2
          },
          {
            label: 'parameters.exposed',
            name: 'e0',
            default: 10000
          },
          {
            label: 'parameters.recovered',
            name: 'r0',
            default: 1000
          }
        ]
      })
    );
  }, 1000);
};

export default slice.reducer;
