import React, {Component} from 'react';
import {connect} from 'react-redux';
import {withTranslation} from 'react-i18next';

import {Spinner} from 'reactstrap';

import GridLayout from './components/GridLayout';

import {initializeApp, setSelected} from './redux/app';
import {setTimeBounds} from './redux/time';

import {RKIDatastore as rki, Tables} from './common/rki-datastore';
import {PopulationDatastore as population} from './common/population-datastore';

import './App.scss';

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
