import React, { Component } from 'react';
import { connect } from 'react-redux';
import { withTranslation } from 'react-i18next';
import Simulation from '../Simulation';
import SEIRChart from '../Graphs/SEIRChart';

import { getSelectedData } from '../../redux/app';
import { getActiveMeasures } from '../../redux/measures';

class Results extends Component {
  render() {
    const { t } = this.props;
    return (
      <>
        <div className="header">{t('results')}</div>
        <div className="content">
          <Simulation />
          <SEIRChart
            {...this.props.seir}
            rki={this.props.rki}
            measures={this.props.measures}
          />
        </div>
      </>
    );
  }
}

const mapState = (state) => {
  return {
    seir: state.seir,
    rki: getSelectedData(state),
    measures: getActiveMeasures(state.measures)
  };
};

export default connect(mapState, {})(withTranslation()(Results));
