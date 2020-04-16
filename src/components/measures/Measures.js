import React, { Component } from 'react';
import { connect } from 'react-redux';
import { withTranslation } from 'react-i18next';
import { CustomInput } from 'reactstrap';




import { addInterval } from '../../redux/measures';

import './Measures.scss';

class Measure extends Component {
    constructor(props) {
        super(props);
        this.state = {
            active: (this.props.measure?.intervals?.length > 0) || false
        };
    }

    toggle() {
        const state = this.state;
        this.setState({
            ...state,
            active: !state.active
        })
    }

    render() {
        const m = this.props.measure;
        const t = this.props.t;
        return (
            <div className="measure px-3">
                <CustomInput 
                    type="switch"
                    id={m.label}
                    value={this.state.active}
                    onChange={this.toggle.bind(this)}
                    className="d-inline"
                />
                <span className="h6 ml-0">{t(m.label)}</span>
            </div>
        );
    }
}

const TranslatedMeasure = withTranslation()(Measure);

class Measures extends Component {
    constructor(props) {
        super(props);
    }

    render() {
        const { t } = this.props;
        return (
            <div className="measures">
                <h4 className="p-2 border-bottom border-secondary">{t('measures.title')}</h4>
                {this.props.measures.map((m, i) => {
                    return <TranslatedMeasure key={i} measure={m} />
                })}
            </div>
        )
    }
}

const mapState = (state) => {
    return {
        measures: state.measures
    };
}

export default connect(mapState, { addInterval })(withTranslation()(Measures))