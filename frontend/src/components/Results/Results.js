import React, {Component} from 'react';
import {connect} from 'react-redux';
import {withTranslation} from 'react-i18next';
import InfectionChart from '../Graphs/InfectionChart';

import {getActiveMeasures} from '../../redux/measures';
import {RKIDatastore as rki} from '../../common/rki-datastore';

import * as numeral from 'numeral';
import {CustomInput} from 'reactstrap';

import './Results.scss';

/**
 * Component for visualising simulation results / rki data
 */
class Results extends Component {
  state = {
    rki: null,
    loading: true,
    logChart: false,
  };

  componentDidMount() {
    if (this.props.selected && this.props.selected.dataset === 'germany') {
      rki.getState(this.props.selected.id).then((data) => {
        this.setState({
          loading: false,
          rki: data,
        });
      });
    }
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.selected?.id !== prevProps.selected?.id) {
      this.setState({
        loading: true,
        rki: null,
      });
      if (['states', 'germany'].includes(this.props.selected.dataset)) {
        rki.getState(this.props.selected.id).then((data) => {
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

  /**
   * Conditional render
   */
  conditionalRender() {
    const {t} = this.props;
    if (this.props.rki === null) {
      return <div>Bitte wählen sie ein Bundesland aus!</div>;
    } else if (this.state.loading) {
      return <div>Loading...</div>;
    } else {
      return (
        <>
          <CustomInput
            type="switch"
            id="log-scale-graph-switch"
            label={t('results.log')}
            checked={this.state.logChart}
            onChange={() => this.setState({logChart: !this.state.logChart})}
            className="p-2 log-switch"
          />
          <InfectionChart
            seir={this.props.seir}
            rki={this.props.rki.all}
            measures={this.props.measures}
            logChart={this.state.logChart}
            style={{height: '100%'}}
          />
        </>
      );
    }
  }

  render() {
    const {t} = this.props;
    return (
      <>
        <div className="header">{t('results.title')}</div>
        <div className="content">
          <div className="top-bar p-1">
            <span className="ml-2">
              {t('selection')} &nbsp;
              {this.props.selected ? this.props.selected.label : t('no-selection')}
              ,&nbsp; {t('population')}:&nbsp;
              {this.props.selected && this.props.selected.population
                ? numeral(this.props.selected.population).format('0,0')
                : '---'}
            </span>
          </div>
          <div className="charts p-1" style={{height: 'calc(100% - 40px)'}}>
            {this.conditionalRender()}
          </div>
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
