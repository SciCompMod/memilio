import React, {Suspense} from 'react';
import ReactDOM from 'react-dom';
import {Provider} from 'react-redux';
import {HashRouter as Router} from 'react-router-dom';
import {Spinner} from 'reactstrap';

import * as numeral from 'numeral';
import * as moment from 'moment';
import 'moment/locale/de';

import App from './App';
import DeveloperPage from './pages/Developer';
import store from './redux/store';
import * as serviceWorker from './serviceWorker';

import './i18n';
import './index.css';

if (process.env.REACT_APP_MODE === 'development') {
  console.warn('USING DEVELOPMENT MODE DO NOT RELEASE THIS!');
}

// setup number library
numeral.register('locale', 'de', {
  delimiters: {
    thousands: '.',
    decimal: ',',
  },
});
numeral.locale('en');

moment.locale('en');

ReactDOM.render(
  <Provider store={store}>
    <Suspense fallback={<Spinner color="white" className="loading-spinner" />}>
      <Router>{process.env.REACT_APP_MODE === 'development' ? <DeveloperPage /> : <App />}</Router>
    </Suspense>
  </Provider>,
  document.getElementById('root')
);

// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: https://bit.ly/CRA-PWA
serviceWorker.unregister();
