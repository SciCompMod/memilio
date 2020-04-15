import { configureStore } from '@reduxjs/toolkit';
import measureReducer from './measures';

export default configureStore({
  reducer: {
    measures: measureReducer,
  },
});
