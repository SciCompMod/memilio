import React, {Component} from 'react';
import {connect} from 'react-redux';
import {withTranslation} from 'react-i18next';

import './AgeChart.scss';

class AgeChart extends Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const {t} = this.props;
    return <></>;
  }
}

const mapState = (state) => {
  return state;
};

export default connect(mapState, null)(withTranslation()(AgeChart));
