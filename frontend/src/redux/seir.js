import { createSlice } from '@reduxjs/toolkit';
import { setSelected } from './app';

const slice = createSlice({
  name: 'seir',
  initialState: {
    startDate: Date.parse('2020-02-24T00:00:00'),
    data: null
  },
  reducers: {
    setStartDate(state, action) {
      state.startDate = action.payload;
    },
    setData(state, action) {
      return {
        ...state,
        data: action.payload
      };
    }
  },
  extraReducers: {
    [setSelected]: (state, action) => {
      state.data = null;
    }
  }
});

export const { setStartDate, setData } = slice.actions;

export default slice.reducer;
