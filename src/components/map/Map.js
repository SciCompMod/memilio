import React, { Component } from 'react';
import { connect } from 'react-redux';

import './Map.scss';

class Map extends Component {

    constructor(props){
        super(props);
        this.node = React.createRef();
    }

    componentDidMount() {
        console.log(this.node.current);
        console.log(this.node.current.clientWidth, this.node.current.clientHeight);
    }

    render() {
        return (
            <div className="map-container">
                <svg ref={this.node}></svg>
            </div>
        )
    }
}

export default connect(null, null)(Map)