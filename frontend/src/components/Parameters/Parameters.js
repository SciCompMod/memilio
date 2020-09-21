import React, {Component} from 'react';
import {withTranslation} from 'react-i18next';
import Parameter from './Parameter';

import './Parameters.scss';

class Parameters extends Component {
  handleChange(update) {
    this.props.onChange(update);
  }

  render() {
    return (
      <div className="parameters mx-3">
        {this.props.parameters.map((p) => {
          return <Parameter key={p.name} parameter={p} onChange={this.handleChange.bind(this)} />;
        })}
      </div>
    );
  }
}

const ParametersTranslated = withTranslation()(Parameters);

export {ParametersTranslated as Parameters};
