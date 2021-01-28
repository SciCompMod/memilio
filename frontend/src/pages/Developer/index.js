import React, {Component} from 'react';

import './styles.scss';
import MainPage from '../Main';
import {Alert} from 'reactstrap';

export default class DeveloperPage extends Component {
  state = {
    visible: true,
  };

  async componentDidMount() {}

  onDismiss() {
    this.setState({
      visible: false,
    });
  }

  render() {
    return (
      <div className="develop-wrapper">
        <Alert
          color="danger"
          isOpen={this.state.visible}
          toggle={this.onDismiss.bind(this)}
          style={{textAlign: 'center'}}
        >
          <b>THIS PAGE IS A DEVELOPMENT VERSION DO NOT RELEASE THIS!</b>
        </Alert>
        <MainPage />
      </div>
    );
  }
}
