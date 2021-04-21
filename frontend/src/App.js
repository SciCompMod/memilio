import React, {Component} from 'react';
import {connect} from 'react-redux';
import {withTranslation} from 'react-i18next';
import am4lang_de_DE from '@amcharts/amcharts4/lang/de_DE';
import * as am4core from '@amcharts/amcharts4/core';
import {Switch, Route, Link} from 'react-router-dom';
import {withRouter} from 'react-router';
import {Spinner} from 'reactstrap';

import {fixUrl} from '~/common/utils';

import Main from './pages/Main';
import About from './pages/About';
import Impressum from './pages/Impressum';
import Dsgvo from './pages/Dsgvo';
import Accessibility from './pages/Accessibility';
import Attribution from './pages/Attribution';

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

  componentWillUnmount() {}

  componentDidUpdate() {
    /* Scroll to top of the page on page update.
       This is necessary when switching between pages like 'impressum' and such.
       
       NOTE: This should only be a temporary solution, since it will probally have a negative impact on
       more advanced components.
    */

    window.scrollTo(0, 0);
  }

  render() {
    const {path, url} = fixUrl(this.props.match);

    return (
      <div className="app">
        <header>
          <div className="title">SARS-CoV-2 Reproduktionszahlen innerhalb Deutschlands</div>
          <div className="logos">
            <div className="hzi-logo">
              <a target="_blank" rel="noopener noreferrer" href="https://www.helmholtz-hzi.de/">
                <img src="assets/logo-hzi2-de.svg" alt="HZI" />
              </a>
            </div>
            <div className="dlr-signet">
              <a target="_blank" rel="noopener noreferrer" href="https://www.dlr.de/">
                <img src="assets/dlr-signet.png" alt="DLR Signet" />
              </a>
            </div>
          </div>
        </header>
        <div className="body">
          <Switch>
            <Route path={`${path}/informationen`}>
              <About />
            </Route>
            <Route path={`${path}/impressum`}>
              <Impressum />
            </Route>
            <Route path={`${path}/datenschutz`}>
              <Dsgvo />
            </Route>
            <Route path={`${path}/barrierefreiheit`}>
              <Accessibility />
            </Route>
            <Route path={`${path}/attribution`}>
              <Attribution />
            </Route>
            <Route path={`${path}/`}>
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
          <Link tabIndex="2" title="Impressum" to={`${url}/impressum`}>
            Impressum
          </Link>
          <Link tabIndex="3" title="Datenschutzerklärung" to={`${url}/datenschutz`}>
            Datenschutzerklärung
          </Link>
          <Link
            tabIndex="4"
            title="Erklärung zur Barrierefreiheit"
            alt="Accessibility statement"
            to={`${url}/barrierefreiheit`}
          >
            Barrierefreiheit
          </Link>
          <Link tabIndex="5" title="Attribution" to={`${url}/attribution`}>
            Attribution
          </Link>
        </div>
      </div>
    );
  }
}

export default connect(null, {})(withTranslation()(withRouter(App)));
