import React, { Component } from 'react';
import { connect } from 'react-redux';
import { Responsive, WidthProvider } from 'react-grid-layout';
import Simulation from './Simulation';
import Map from './Map';
import Parameters from './Parameters';
import Measures from './Measures';
import SEIRChart from './Graphs/SEIRChart';

import { getActiveMeasures } from '../redux/measures';

import './GridLayout.scss';

const ResponsiveGridLayout = WidthProvider(Responsive);

class ResponsiveGrid extends Component {
  render() {
    // {lg: layout1, md: layout2, ...}
    //const layouts = getLayoutsFromSomewhere();
    return (
      <>
        <div className="earth">
          <img src="assets/earth.png" alt="Welt" />
        </div>
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
        <ResponsiveGridLayout
          className="layout position-relative"
          breakpoints={{ lg: 1200, md: 996, sm: 768, xs: 480, xxs: 0 }}
          cols={{ lg: 12, md: 10, sm: 6, xs: 4, xxs: 2 }}
          rowHeight={100}
          containerPadding={[10, 10]}
          margin={[5, 5]}
          draggableHandle=".grid-draggable-handle"
        >
          <div
            key="map"
            data-grid={{ x: 0, y: 0, w: 5, h: 4 }}
            className="grid-box"
          >
            <i className="fas fa-arrows-alt grid-draggable-handle"></i>
            <Map />
          </div>
          <div
            key="measures"
            data-grid={{ x: 0, y: 4, w: 3, h: 4 }}
            className="grid-box"
          >
            <i className="fas fa-arrows-alt grid-draggable-handle"></i>
            <Measures />
          </div>
          <div
            key="parameters"
            data-grid={{ x: 3, y: 4, w: 2, h: 4 }}
            className="grid-box"
          >
            <i className="fas fa-arrows-alt grid-draggable-handle"></i>
            <Parameters />
            <Simulation />
          </div>
          <div
            key="seir-chart"
            data-grid={{ x: 5, y: 0, w: 7, h: 5 }}
            className="grid-box"
          >
            <i className="fas fa-arrows-alt grid-draggable-handle"></i>
            <SEIRChart {...this.props.seir} measures={this.props.measures} />
          </div>
        </ResponsiveGridLayout>
      </>
    );
  }
}

const mapState = (state) => {
  return {
    seir: state.seir,
    measures: getActiveMeasures(state.measures)
  };
};

export default connect(mapState, {})(ResponsiveGrid);
