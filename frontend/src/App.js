import React, {Component} from 'react';
import {connect} from 'react-redux';
import {withTranslation} from 'react-i18next';
import am4lang_de_DE from '@amcharts/amcharts4/lang/de_DE';
import * as am4core from '@amcharts/amcharts4/core';

import {Spinner} from 'reactstrap';

import GridLayout from './components/GridLayout';

import {initializeApp, setSelected} from './redux/app';
import {setTimeBounds} from './redux/time';

import {RKIDatastore as rki, Tables} from './common/rki-datastore';
import {PopulationDatastore as population} from './common/population-datastore';
import './common/rki-sql-store';

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
    };
  }

  async componentDidMount() {
    await this.props.initializeApp();
    rki.populate([Tables.STATES, Tables.COUNTIES]).then(() => {
      Promise.all([rki.getTimeBounds(Tables.STATES), population.getByKey(-1)])
        .then((r) => {
          console.log(r);
          const [{start, end}, p] = r;
          this.props.setTimeBounds({
            start,
            end,
          });
          console.log(p);
          this.props.setSelected({
            id: 0,
            dataset: 'germany',
            label: 'Deutschland',
            population: p.population,
          });
          this.setState({loading: false});
        })
        .catch((err) => console.log(err));
    });
  }

  render() {
    if (this.state.loading) {
      return <Spinner color="white" className="loading-spinner" />;
    }

    return <GridLayout />;
  }
}

export default connect(null, {initializeApp, setTimeBounds, setSelected})(withTranslation()(App));
