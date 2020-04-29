import React, { Component } from 'react';
import { connect } from 'react-redux';
import { withTranslation } from 'react-i18next';
import { updateParameter } from '../../redux/parameters';
import { Label, Input, CustomInput, FormGroup, Col } from 'reactstrap';

import './Parameters.scss';

class Parameter extends Component {
  constructor(props) {
    super(props);
    this.timeout = null;
    this.state = {
      value: this.props.parameter.default
    };
  }

  componentWillUnmount() {
    if (this.timeout) {
      clearTimeout(this.timeout);
      this.timeout = null;
    }
  }

  update(value) {
    this.setState({
      value
    });

    if (this.timeout) {
      clearTimeout(this.timeout);
    }

    this.timeout = setTimeout(() => {
      this.timeout = null;

      this.props.updateParameter({
        ...this.props.parameter,
        value: this.state.value
      });
    }, 300);
  }

  renderInput() {
    const { parameter } = this.props;
    switch (parameter.type) {
      case 'slider':
        return (
          <FormGroup row className="m-0 mb-2">
            <Col xs="2" className="p-0 pr-1">
              <div className="border border-dark rounded-lg p-1 text-center">
                {this.state.value}
              </div>
            </Col>
            <Col xs="10" className="p-0 pl-1">
              <CustomInput
                type="range"
                value={this.state.value}
                id={parameter.name}
                name={parameter.name}
                min={parameter.min || 0}
                max={parameter.max || 10}
                step={parameter.step || 1}
                onChange={(event) =>
                  this.update(parseFloat(event.target.value))
                }
                tooltip="off"
                size="sm"
              />
            </Col>
          </FormGroup>
        );
      default:
        return (
          <FormGroup className="m-0 mb-2">
            <Input
              value={this.state.value}
              name={parameter.name}
              onChange={(event) => this.update(parseFloat(event.target.value))}
            />
          </FormGroup>
        );
    }
  }

  render() {
    const { t } = this.props;
    return (
      <div className="parameter">
        <Label className="h6 m-0">{t(this.props.parameter.label)}</Label>
        {this.renderInput()}
      </div>
    );
  }
}

const TransletedParameter = connect(null, { updateParameter })(
  withTranslation()(Parameter)
);

class Parameters extends Component {
  constructor(props) {
    super(props);
    this.state = {
      parameters: []
    };
  }

  static getDerivedStateFromProps(nextProps, prevState) {
    const { parameters } = nextProps;
    return {
      ...prevState,
      parameters
    };
  }

  handleChange(event) {
    console.log(event.target.name, event.target.value);
  }

  handleSubmit(event) {
    event.preventDefault();
  }

  render() {
    const { t } = this.props;
    return (
      <div className="parameters">
        <h5 className="p-2 border-bottom border-secondary">
          {t('parameters.title')}
        </h5>
        <form className="parameter" onSubmit={this.handleSubmit}>
          {this.state.parameters.map((p) => {
            return <TransletedParameter key={p.name} parameter={p} />;
          })}
        </form>
      </div>
    );
  }
}

const mapState = (state) => {
  return {
    parameters: state.parameters.parameters
  };
};

export default connect(mapState, null)(withTranslation()(Parameters));
