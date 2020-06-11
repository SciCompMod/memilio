import { createSlice } from '@reduxjs/toolkit';

const slice = createSlice({
    name: 'time',
    initialState: {
        time: Date.parse('2020-02-24T00:00:00'),
    },
    reducers: {
        setCurrentTime(state, action) {
            state.time = action.payload;
        }
    }
});

export const { setCurrentTime } = slice.actions;

export default slice.reducer;
