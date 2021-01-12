import React, {Component} from 'react';

import './styles.scss';

export default class DevoloperPage extends Component {
  state = {};

  async componentDidMount() {}

  render() {
    return (
      <div className="develop-wrapper">
        THis is a wrapper
        {this.props.children}
      </div>
    );
  }
}
