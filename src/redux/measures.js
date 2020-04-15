import { createSlice } from '@reduxjs/toolkit';

const slice = createSlice({
    name: 'measures',
    initialState: [],
    reducers: {
        addInterval: (state, action) => {
            console.log(state, action);
            const measure = state.find(m => m.id === action.payload.measure);
            if (measure) {
                measure.intervals.push(action.payload.interval);
            }
        }
    },
});

export const { addInterval } = slice.actions;

export default slice.reducer;
