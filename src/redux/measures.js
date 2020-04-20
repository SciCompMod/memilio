import { createSlice } from '@reduxjs/toolkit';
import { init } from './app';
import { v4 as uuid } from 'uuid';

const findById = (list, id) => {
    return list.find(m => m.id === id);
}

const sortByDate = (list) => {
    return list.sort((a, b) => {
        if (a.start === b.start) {
            return a.end - b.end;
        }
        return a.start - b.start;
    });
}

const slice = createSlice({
    name: 'measures',
    initialState: [],
    reducers: {
        addInterval: (state, action) => {
            // search measure
            const measure = findById(state, action.payload.measure);
            if (measure) {
                const id = uuid();
                const {start, end} = action.payload;
                measure.intervals.push({ 
                    id,
                    start, 
                    end 
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
                if(interval) {
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
                    return Object.assign({
                        id: uuid(),
                        intervals: [],
                        damping: 0,
                        label: ""
                    }, m);
                });
            }
            return state;
        }
    }
});

export const { addInterval, editInterval } = slice.actions;

export default slice.reducer;
