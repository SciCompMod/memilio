import { createSlice } from '@reduxjs/toolkit';
import { init } from './app';
import { v4 as uuid } from 'uuid';

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

const slice = createSlice({
  name: 'measures',
  initialState: [],
  reducers: {
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
    addInterval: (state, action) => {
      // search measure
      const measure = findById(state, action.payload.measure);
      if (measure) {
        const id = uuid();
        const { start, end } = action.payload;
        measure.intervals.push({
          id,
          start,
          end,
          active: true
        });

        // reorder by start date and end date
        measure.intervals = sortByDate(measure.intervals);
      }
    },
    editInterval: (state, action) => {
      // search measure
      const measure = findById(state, action.payload.measure);
      if (measure) {
        // search interval
        const interval = findById(measure.intervals, action.payload.id);
        if (interval) {
          interval.start = action.payload.start;
          interval.end = action.payload.end;
        }

        // reorder by start date and end date
        measure.intervals = sortByDate(measure.intervals);
      }
    }
  },
  extraReducers: {
    [init]: (state, action) => {
      if (action.payload.measures) {
        return action.payload.measures.map((m, i) => {
          if (m.intervals) {
            m.intervals = m.intervals.map((i) => ({
              ...i,
              id: uuid(),
              active: true
            }));
            m.active = true;
          }
          return Object.assign(
            {
              id: uuid(),
              intervals: [],
              damping: 0,
              label: '',
              active: false
            },
            m
          );
        });
      }
      return state;
    }
  }
});

export const {
  addInterval,
  editInterval,
  activateMeasure,
  deactivateMeasure
} = slice.actions;

export const getActiveMeasures = (state) => {
  return state
    .map((m) => {
      if (!m.active) {
        return null;
      }

      return {
        ...m,
        intervals: m.intervals
          .map((i) => (i.active ? i : null))
          .filter((i) => i !== null)
      };
    })
    .filter((m) => m !== null);
};

export default slice.reducer;
