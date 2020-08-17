import {createSlice} from '@reduxjs/toolkit';
import {roundToUTCMidnight} from '../common/utils';

const START_DATE = Date.parse('2020-02-24T00:00:00');

const slice = createSlice({
  name: 'time',
  initialState: {
    startDate: START_DATE,
    endDate: START_DATE,
    currentDate: START_DATE,
  },
  reducers: {
    setStartDate(state, action) {
      state.startDate = action.payload;
    },
    setCurrentDate(state, action) {
      state.currentDate = roundToUTCMidnight(action.payload);
    },
    setEndDate(state, action) {
      state.endDate = action.payload;
    },
    setTimeBounds(state, action) {
      state.startDate = action.payload.start;
      state.currentDate = action.payload.start;
      state.endDate = action.payload.end;
    },
  },
});

export const {setCurrentDate, setStartDate, setEndDate, setTimeBounds} = slice.actions;

export default slice.reducer;
