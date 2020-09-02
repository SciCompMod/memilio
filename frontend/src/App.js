import React, {Component} from 'react';
import {connect} from 'react-redux';
import {withTranslation} from 'react-i18next';
import am4lang_de_DE from '@amcharts/amcharts4/lang/de_DE';
import * as am4core from '@amcharts/amcharts4/core';

import GridLayout from './components/GridLayout';

import {initializeApp, fetchData} from './redux/app';
import './App.scss';

// AmCharts defaults to English as a locale and not the browser default,
// so we have to set it manually.
if (navigator.language.includes('de')) {
  am4core.options.defaultLocale = am4lang_de_DE;
}

class App extends Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  componentDidMount() {
    this.props.initializeApp();
    this.props.fetchData();
  }

  render() {
    return <GridLayout />;
  }

  /* render() {
    const { t } = this.props;
    return (
      <Container fluid className="wrapper">
        <header>
          <span className="logo">
            <img
              src="assets/logo.png"
              alt="Deutsches Zentrum für Luft- und Raumfahrt (DLR)"
            />
            <span className="logo-text">
              Deutsches Zentrum
              <br />
              für Luft- und Raumfahrt
            </span>
          </span>
        </header>
        <Row className="main-panel">
          <Col xl="5" lg="5" md="12" sm="12" xs="12">
            <Row className="mb-1 h-100">
              <Col className="map-container p-0">
                <Map />
              </Col>
              <div className="w-100 d-lg-block d-md-none"></div>
              <Col className="p-0 settings-container mt-1">
                <Nav tabs className="d-lg-none">
                  <NavItem>
                    <NavLink
                      onClick={() => {
                        this.toggle(1);
                      }}
                    >
                      {t('measures.title')}
                    </NavLink>
                  </NavItem>
                  <NavItem>
                    <NavLink
                      onClick={() => {
                        this.toggle(2);
                      }}
                    >
                      {t('parameters.title')}
                    </NavLink>
                  </NavItem>
                </Nav>
                <TabContent
                  activeTab={`${this.state.active}`}
                  className="row h-100"
                >
                  <TabPane tabId="1" className="col-lg-6 d-lg-block">
                    <Row>
                      <Col className="pr-0">
                        <Measures />
                      </Col>
                    </Row>
                  </TabPane>
                  <TabPane
                    tabId="2"
                    className="col-lg-6 d-lg-block border-secondary border-left"
                  >
                    <Row>
                      <Col className="pl-0">
                        <Parameters />
                      </Col>
                    </Row>
                  </TabPane>
                </TabContent>
              </Col>
            </Row>
          </Col>
          <Col xl="7" lg="7" md="12" sm="12" xs="12" className="px-1">
            <div className="graphs-container bg-primary">
              <Graphs data={this.data} />
            </div>
          </Col>
        </Row>
      </Container>
    );
  } */
}

export default connect(null, {initializeApp, fetchData})(withTranslation()(App));
