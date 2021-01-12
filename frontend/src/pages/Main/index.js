import React, {Component} from 'react';

import './styles.scss';

class MainPage extends Component {
  render() {
    return (
      <div className="main">
        <div className="left">Here goes the map</div>
        <div className="right">
          <div className="graph">Here goes the graph</div>
        </div>
      </div>
    );
  }
}

export default MainPage;
