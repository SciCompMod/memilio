import React, {Component} from 'react';
import {connect} from 'react-redux';
import {withTranslation} from 'react-i18next';
import Simulation from '../Simulation';
import Chart from '../Graphs/Chart';

import {getSelectedData} from '../../redux/app';
import {getActiveMeasures} from '../../redux/measures';

import * as numeral from 'numeral';
import {CustomInput} from 'reactstrap';

class Results extends Component {
  state = {
    logChart: false,
  };

  _render() {
    if (this.props.rki === null) {
      return <div>Bitte wählen sie ein Bundesland aus!</div>;
    }
    return (
      <Chart
        seir={this.props.seir}
        rki={this.props.rki.all}
        measures={this.props.measures}
        logChart={this.state.logChart}
        style={{height: '100%'}}
      />
    );
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
            <span className="h3 ml-0"> | Logarithmic: </span>
            <CustomInput
              type="switch"
              id="log-scale-graph-switch"
              checked={this.state.logChart}
              onChange={() => this.setState({logChart: !this.state.logChart})}
              className="d-inline"
            />
          </div>
          <div className="charts p-1" style={{height: 'calc(100% - 39px)'}}>
            {this._render()}
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
    rki: getSelectedData(state),
    measures: getActiveMeasures(state.measures),
  };
};

const ResultsTranslated = withTranslation()(Results);
const ResultsConnected = connect(mapState, {})(ResultsTranslated);

export {ResultsConnected as Results};
