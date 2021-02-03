import React, {Component} from 'react';
import {Button, Alert} from 'reactstrap';
import {Link, Redirect, Route, Switch} from 'react-router-dom';

import App from '../../App';

import './styles.scss';

export default class DeveloperPage extends Component {
  state = {
    open: false,
    small: true,
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
          className="developer-alert"
          color="danger"
          isOpen={this.state.visible}
          toggle={this.onDismiss.bind(this)}
          style={{textAlign: 'center'}}
        >
          <b>This is a development version. DO NOT RELEASE THIS!</b>
        </Alert>
        <div className={`sidebar ${this.state.small ? 'small' : 'normal'} ${this.state.open ? 'open' : 'closed'}`}>
          <header>
            <Button
              color="link"
              onClick={() => {
                this.setState({
                  small: !this.state.small,
                });
              }}
            >
              <i className="fas fa-angle-double-left"></i>
              <i className="fas fa-angle-double-right"></i>
            </Button>
          </header>
          <ul
            className="menu"
            onMouseEnter={() => {
              this.setState({
                open: true,
              });
            }}
            onMouseLeave={() => {
              this.setState({
                open: false,
              });
            }}
          >
            <li className="menu-item active">
              <Link to="/production">
                <i className="fa fa-chart-pie"></i>Production
              </Link>
            </li>
            <li className="menu-item">
              <span className="more">more links here...</span>
            </li>
          </ul>
        </div>
        <div className="body">
          <main>
            <Switch>
              <Route path="/production">
                <App />
              </Route>
              <Redirect from="/" to="/production" />
            </Switch>
          </main>
        </div>
      </div>
    );
  }
}
