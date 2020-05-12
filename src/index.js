import React, { Suspense } from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';
import App from './App';
import store from './redux/store';
import * as serviceWorker from './serviceWorker';

import { Spinner } from 'reactstrap';

import './i18n';

import './index.css';

ReactDOM.render(
  <Provider store={store}>
    <Suspense fallback={<Spinner color="white" className="loading-spinner" />}>
      <App />
    </Suspense>
  </Provider>,
  document.getElementById('root')
);

// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: https://bit.ly/CRA-PWA
serviceWorker.unregister();
