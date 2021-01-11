import React, {Component} from 'react';
import {connect} from 'react-redux';
import {withTranslation} from 'react-i18next';
import am4lang_de_DE from '@amcharts/amcharts4/lang/de_DE';
import * as am4core from '@amcharts/amcharts4/core';
import {BrowserRouter as Router, Switch, Route} from 'react-router-dom';

import {Spinner} from 'reactstrap';

import Main from './pages/Main';

import {initializeApp, setSelected} from './redux/app';
import {setTimeBounds} from './redux/time';

import rki from './common/datastore/sql/rki-sql-store';
import population from './common/datastore/idb/population-datastore';

import './App.scss';

// AmCharts defaults to English as a locale and not the browser default,
// so we have to set it manually.
if (navigator.language.includes('de')) {
  am4core.options.defaultLocale = am4lang_de_DE;
}

class App extends Component {
  constructor(props) {
    super(props);
    this.state = {
      loading: true,
      label: '',
    };
  }

  progress(state) {
    return new Promise((resolve, reject) => {
      this.setState(state, () => {
        resolve();
      });
    });
  }

  delay(timeout) {
    return new Promise((resolve, reject) => {
      setTimeout(() => {
        resolve();
      }, timeout);
    });
  }

  async componentDidMount() {
    await this.progress({label: 'Loading configuration'});
    await this.props.initializeApp();

    await this.progress({label: 'Loading data'});
    await population.populate();
    await rki.populate();

    await this.progress({label: 'Initializing application'});
    const p = await population.getByKey('0');
    const {start, end} = await rki.getTimeBounds();

    this.props.setTimeBounds({start, end});
    this.props.setSelected({
      id: '0',
      dataset: 'germany',
      label: 'Deutschland',
      population: p.population,
    });

    await this.progress({label: 'Done'});
    this.setState({loading: false});
  }

  render() {
    return (
      <Router>
        <Switch>
          <Route path="/">
            <Main />
          </Route>
        </Switch>
        {this.state.loading ? (
          <div className="overlay">
            <div className="loading-progress">
              <div className="label">{this.state.label}</div>
              <Spinner color="primary" />
            </div>
          </div>
        ) : (
          <></>
        )}
      </Router>
    );
  }
}

export default connect(null, {initializeApp, setTimeBounds, setSelected})(withTranslation()(App));
