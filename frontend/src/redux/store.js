import {configureStore} from '@reduxjs/toolkit';
import measureReducer from './measures';
import appReducer from './app';
import parameterReducer from './parameters';
import modelsReducer from './models';
import seirReducer from './seir';
import timeReducer from './time';

import thunk from 'redux-thunk';

export default configureStore({
  reducer: {
    measures: measureReducer,
    parameters: parameterReducer,
    models: modelsReducer,
    seir: seirReducer,
    time: timeReducer,
    app: appReducer,
  },
  middleware: [thunk],
});
