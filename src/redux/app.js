import { createSlice } from '@reduxjs/toolkit';

const slice = createSlice({
    name: 'app',
    initialState: {},
    reducers: {
        init: (state, action) => {
            
        }
    },
});

export const { init } = slice.actions;

export const initialize = () => dispatch => {
    setTimeout(() => {
        dispatch(init({
            measures: [{
                id: 1,
                intervals: []
            }],
            parameters: {
                a: 5
            }
        }));
    }, 1000);
};

export default slice.reducer;
