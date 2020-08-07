import React, {Component} from 'react';
import {connect} from 'react-redux';

import InteractiveMap from '../../common/interactive_map';
import {setSelected} from '../../redux/app';

import './Map.scss';

class GermanyMap extends Component {
  constructor(props) {
    super(props);
    this.node = React.createRef();
  }

  componentDidMount() {
    this.map = new InteractiveMap(this.node.current);
    this.map.onSelect((selected) => {
      const {dataset, id} = selected;
      if (!dataset || !id) {
        this.props.setSelected(null);
      } else {
        this.props.setSelected(selected);
      }
    });
  }

  render() {
    return (
      <div className="map">
        <svg ref={this.node}></svg>
      </div>
    );
  }
}

const ConnectedMap = connect(null, {setSelected})(GermanyMap);

export {ConnectedMap as GermanyMap};
