import { createSlice } from '@reduxjs/toolkit';

const slice = createSlice({
  name: 'seir',
  initialState: {
    startDate: Date.parse('2020-02-24'),
    S: [],
    E: [],
    I: [],
    R: []
  },
  reducers: {
    setStartDate(state, action) {
      state.startDate = action.payload;
    },
    setData(state, action) {
      const { S, E, I, R } = action.payload;
      return {
        ...state,
        S,
        E,
        I,
        R
      };
    }
  },
  extraReducers: {}
});

export const { setStartDate, setData } = slice.actions;

export default slice.reducer;
