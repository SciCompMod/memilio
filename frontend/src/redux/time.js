import { createSlice } from '@reduxjs/toolkit';

const START_DATE = Date.parse('2020-02-24T00:00:00');

const slice = createSlice({
    name: 'time',
    initialState: {
        startDate: START_DATE,
        endDate: START_DATE,
        currentDate: START_DATE
    },
    reducers: {
        setStartDate(state, action) {
            state.startDate = action.payload;
        },
        setCurrentDate(state, action) {
            state.currentDate = action.payload;
        },
        setEndDate(state, action) {
            state.endDate = action.payload
        }
    }
});

export const { setCurrentDate, setStartDate, setEndDate } = slice.actions;

export default slice.reducer;
