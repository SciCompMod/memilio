import React, { Component } from 'react';
import { connect } from 'react-redux';
import { addInterval } from '../redux/measures';

class Measures extends Component {
    render() {
        return (
            <div>
                <i
                    onClick={() => this.props.addInterval({
                        measure: 1,
                        interval: {
                            start: 0,
                            end: 0
                        }
                    })}
                    className="fa fa-plus"
                ></i>
            </div>
        )
    }
}

const mapState = (state) => {
    console.log(state);
    return {
        measure: 2
    };
}

export default connect(mapState, { addInterval })(Measures)