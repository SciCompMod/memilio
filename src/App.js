import React, { Component } from 'react';
import { Container, Row, Col } from 'reactstrap';
import Map from './components/map/Map'

import './App.scss';

class App extends Component {
  render() {
    return (
      <Container fluid className="wrapper">
        <Row className="main-panel">
          <Col xs="4" className="bg-primary">
            <Map />
          </Col>
          <Col xs="6" className="bg-secondary"></Col>
          <Col xs="2" className="bg-success"></Col>
        </Row>
      </Container>
    );
  }
}

export default App;
