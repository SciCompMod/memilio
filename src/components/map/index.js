import React, { Component } from 'react';
import { connect } from 'react-redux';

import InteractiveMap from '../../common/interactive_map';

import './Map.scss';

class Map extends Component {
  constructor(props) {
    super(props);
    this.node = React.createRef();
  }

  componentDidMount() {
    this.map = new InteractiveMap(this.node.current);
  }

  render() {
    return (
      <div className="map">
        <svg ref={this.node}></svg>
      </div>
    );
  }
}

export default connect(null, null)(Map);
