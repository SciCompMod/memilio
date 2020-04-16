import { createSlice } from '@reduxjs/toolkit';
import {init} from './app';

function clone(object) {
    return Object.assign({}, object);
}

const slice = createSlice({
    name: 'parameters',
    initialState: {
        defaults: [],
        parameters: []
    },
    reducers: {},
    extraReducers: {
        [init]: (state, action) => {
            if(action.payload.parameters) {
                const defaults = clone({defaults: action.payload.parameters});
                const parameters = clone({parameters: action.payload.parameters});
                return Object.assign({}, state, defaults, parameters);
            }
            return state;
        }
    }
});


export default slice.reducer;
