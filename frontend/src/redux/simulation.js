/*
    This represents the currently displayed simulation data.
    It includes:
        - the used simulation model
        - all set parameters
        - all configured measures
*/

import {createSlice} from '@reduxjs/toolkit';
import {init} from './app';
import {v4 as uuid} from 'uuid';
import * as dayjs from 'dayjs';

const findById = (list, id) => {
  return list.find((m) => m.id === id);
};

const sortByDate = (list) => {
  return list.sort((a, b) => {
    if (a.start === b.start) {
      return a.end - b.end;
    }
    return a.start - b.start;
  });
};

/**
 * Utility function to get milliseconds from beginning of current date
 */
const today = () => {
  return dayjs().startOf('day').valueOf();
};

/**
 * Utility function to get milliseconds for last day of year
 */
const endOfYear = () => {
  return dayjs().endOf('year').startOf('day').valueOf();
};

/**
 * Sanitize start and end date.
 *
 * @param {start, end} interval
 */
const sanitizeInterval = (interval) => {
  // check start date
  if (!interval.start) {
    // if start date is undefined, null or 0, use current day
    interval.start = today();
  } else if (typeof interval.start === 'string') {
    // if start date is a string, convert it to milliseconds
    interval.start = Date.parse(interval.start);
  } else if (typeof interval.start === 'number') {
    // if start date is a number, assume it already is in milliseconds
    // do nothing
  }

  // check end date
  if (!interval.end) {
    // if end date is undefined, null or 0, use end of year
    interval.end = endOfYear();
  } else if (typeof interval.end === 'string') {
    // if end date is a string, convert it to milliseconds
    interval.end = Date.parse(interval.end);
  } else if (typeof interval.end === 'number') {
    // if end date is a number, assume it already is in milliseconds
    // do nothing
  }
  return interval;
};

const slice = createSlice({
  name: 'simulation',
  initialState: {
    metadata: null, // meta data, such as which state / county / population ...
    start: null, // start date of simulation
    days: 0, // length of simulation
    model: null, // simulation model used
    measures: [], // set measures
    parameters: [], // set parameters
    results: null, // simulation result (TODO: store in IndexedDB or sql.js)
  },
  reducers: {
    addMeasure: (state, action) => {
      const measure = sanitizeInterval(action.payload);

      state.measures.push({
        ...measure,
        id: uuid,
        active: true,
      });

      state.measures = sortByDate(state.measures);
    },
    editMeasure: (state, action) => {
      // search measure
      const measure = findById(state, action.payload.id);
      if (measure) {
        const sanitized = sanitizeInterval(action.payload);

        measure.type = sanitized.type;
        measure.start = sanitized.start;
        measure.end = sanitized.end;
        measure.active = sanitized.active;

        // reorder by start date and end date
        state.measures = sortByDate(state.measures);
      }
    },
    activateMeasure: (state, action) => {
      const measure = findById(state, action.payload.id);
      if (measure) {
        measure.active = true;
      }
    },
    deactivateMeasure: (state, action) => {
      const measure = findById(state, action.payload.id);
      if (measure) {
        measure.active = false;
      }
    },
  },
  extraReducers: {
    [init]: (state, action) => {
      return state;
    },
  },
});

export const {addMeasure, editMeasure, activateMeasure, deactivateMeasure} = slice.actions;

export const getActiveMeasures = (state) => {
  return state.measures.reduce((acc, measure) => {
    if (measure.active) {
      acc.push(measure);
    }
    return acc;
  }, []);
};

export default slice.reducer;
