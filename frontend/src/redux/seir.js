import { createSlice } from '@reduxjs/toolkit';
import { setSelected } from './app';

const slice = createSlice({
  name: 'seir',
  initialState: {
    startDate: Date.parse('2020-02-24T00:00:00'),
    data: null,
    regionData: null
  },
  reducers: {
    setStartDate(state, action) {
      state.startDate = action.payload;
    },
    setData(state, action) {
      state.data = action.payload;
    },
    setRegionData(state, action) {
      state.regionData = action.payload;
    }
  },
  extraReducers: {
    [setSelected]: (state, action) => {
      state.data = null;
      state.regionData = null;
    }
  }
});

export const { setStartDate, setData, setRegionData } = slice.actions;

export default slice.reducer;
