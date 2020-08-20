import React, {Component} from 'react';
import {connect} from 'react-redux';
import {withTranslation} from 'react-i18next';
import Simulation from '../Simulation';
import SEIRChart from '../Graphs/SEIRChart';

import {getActiveMeasures} from '../../redux/measures';
import {RKIDatastore as rki} from '../../common/rki-datastore';

import * as numeral from 'numeral';

class Results extends Component {
  state = {
    rki: null,
    loading: true,
  };

  componentDidMount() {
    console.log('result did mount', this.props.selected);
    if (this.props.selected && this.props.selected.dataset === 'germany') {
      rki.getState(this.props.selected.id).then((data) => {
        console.log('loaded', data);
        this.setState({
          loading: false,
          rki: data,
        });
      });
    }
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.selected?.id !== prevProps.selected?.id) {
      console.log('result did update');
      console.log(this.props.selected, prevProps.selected);
      this.setState({
        loading: true,
        rki: null,
      });
      if (['states', 'germany'].includes(this.props.selected.dataset)) {
        rki.getState(this.props.selected.id).then((data) => {
          console.log('loaded', data);
          this.setState({
            loading: false,
            rki: data,
          });
        });
      } else if (this.props.selected.dataset === 'counties') {
        rki.getCounty(this.props.selected.id).then((data) => {
          this.setState({
            loading: false,
            rki: data,
          });
        });
      }
    }
  }

  _render() {
    if (!this.props.selected) {
      return <div>Bitte wählen sie ein Bundesland aus!</div>;
    } else if (this.state.loading) {
      return <div>Loading...</div>;
    } else {
      return <SEIRChart seir={this.props.seir} rki={this.state.rki} measures={this.props.measures} />;
    }
  }

  render() {
    const {t} = this.props;
    return (
      <>
        <div className="header">{t('results')}</div>
        <div className="content">
          <div className="top-bar p-1">
            <Simulation />
            <span className="ml-2">
              {t('selection')} &nbsp;
              {this.props.selected ? this.props.selected.label : t('no-selection')}
              ,&nbsp; {t('population')}:&nbsp;
              {this.props.selected && this.props.selected.population
                ? numeral(this.props.selected.population).format('0,0')
                : '---'}
            </span>
          </div>
          <div className="charts p-1">{this._render()}</div>
        </div>
      </>
    );
  }
}

const mapState = (state) => {
  return {
    selected: state.app.selected,
    seir: state.seir.data,
    rki: null, //getSelectedData(state),
    measures: getActiveMeasures(state.measures),
  };
};

const ResultsTranslated = withTranslation()(Results);
const ResultsConnected = connect(mapState, {})(ResultsTranslated);

export {ResultsConnected as Results};
