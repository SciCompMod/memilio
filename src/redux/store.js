import { configureStore } from '@reduxjs/toolkit';
import measureReducer from './measures';
import appReducer from './app';
import parameterReducer from './parameters';

export default configureStore({
  reducer: {
    measures: measureReducer,
    parameters: parameterReducer,
    app: appReducer
  },
});
