import {createSlice} from '@reduxjs/toolkit';
import {init} from './app';

function clone(object) {
  return Object.assign({}, object);
}

const slice = createSlice({
  name: 'parameters',
  initialState: {
    defaults: [],
    parameters: [],
  },
  reducers: {
    updateParameter: (state, action) => {
      const parameter = state.parameters.find((p) => p.name === action.payload.name);
      if (parameter) {
        parameter.value = action.payload.value;
      }
    },
  },
  extraReducers: {
    [init]: (state, action) => {
      if (action.payload.parameters) {
        const defaults = clone({defaults: action.payload.parameters});
        const parameters = clone({parameters: action.payload.parameters});
        return Object.assign({}, state, defaults, parameters);
      }
      return state;
    },
  },
});

export const {updateParameter} = slice.actions;

export const getParameterMap = (state) => {
  let map = {};

  state.parameters.forEach((p) => {
    map[p.name] = p.value || p.default;
  });

  return map;
};

export default slice.reducer;
