import { createSlice } from '@reduxjs/toolkit';
import {init} from './app';

const slice = createSlice({
    name: 'measures',
    initialState: [],
    reducers: {
        addInterval: (state, action) => {
            console.log('addInterval', state, action);
            const measure = state.find(m => m.id === action.payload.measure);
            if (measure) {
                measure.intervals.push(action.payload.interval);
            }
        }
    },
    extraReducers: {
        [init]: (state, action) => {
            if(action.payload.measures) {
                return action.payload.measures;
            }
            return state;
        }
    }
});

export const { addInterval } = slice.actions;

export default slice.reducer;
