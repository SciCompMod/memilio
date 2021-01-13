import React, {Component} from 'react';
import {connect} from 'react-redux';
import {withTranslation} from 'react-i18next';
import am4lang_de_DE from '@amcharts/amcharts4/lang/de_DE';
import * as am4core from '@amcharts/amcharts4/core';
import {HashRouter as Router, Switch, Route, Link} from 'react-router-dom';

import {Spinner} from 'reactstrap';

import Main from './pages/Main';

import './App.scss';
import About from './pages/About';
import Impressum from './pages/Impressum';
import Dsgvo from './pages/Dsgvo';
import Accessibility from './pages/Accessibility';

// AmCharts defaults to English as a locale and not the browser default,
// so we have to set it manually.
if (navigator.language.includes('de')) {
  am4core.options.defaultLocale = am4lang_de_DE;
}

class App extends Component {
  constructor(props) {
    super(props);
    this.state = {
      loading: false,
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

  async componentDidMount() {}

  render() {
    return (
      <div className="app">
        <Router>
          <header>
            <div className="title">SARS-CoV-2 Reproduktionszahlen innerhalb Deutschlands</div>
            <div className="logos">
              <div className="hzi-logo">
                <img src="assets/logo-hzi2-de.svg" alt="HZI" />
              </div>
              <div className="dlr-signet">
                <img src="assets/dlr-signet.png" alt="DLR Signet" />
              </div>
            </div>
          </header>
          <div className="body">
            <Switch>
              <Route path="/impressum">
                <Impressum />
              </Route>
              <Route path="/datenschutz">
                <Dsgvo />
              </Route>
              <Route path="/informationen">
                <About />
              </Route>
              <Route path="/barrierefreiheit">
                <Accessibility />
              </Route>
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
          </div>
          <div className="links">
            <Link tabIndex="2" title="Weitere Informationen zu der Webseite" to="/informationen">
              Informationen
            </Link>
            <Link tabIndex="3" title="Impressum" to="/impressum">
              Impressum
            </Link>
            <Link tabIndex="4" title="Datenschutzerklärung" to="/datenschutz">
              Datenschutzerklärung
            </Link>
            <Link
              tabIndex="5"
              title="Erklärung zur Barrierefreiheit"
              alt="Accessibility statement"
              to="/barrierefreiheit"
            >
              Barrierefreiheit
            </Link>
          </div>
        </Router>
      </div>
    );
  }
}

export default connect(null, {})(withTranslation()(App));
