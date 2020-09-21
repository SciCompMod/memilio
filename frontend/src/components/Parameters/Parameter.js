import React, { Component } from 'react';
import { withTranslation } from 'react-i18next';
import { Label, Input, CustomInput, FormGroup, Col } from 'reactstrap';

class Parameter extends Component {
  constructor(props) {
    super(props);
    this.timeout = null;
    this.state = {
      value: this.props.parameter.value || this.props.parameter.default
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
      this.props.onChange({
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
              <div className="h4 border border-dark rounded-lg p-1 text-center">
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
        <Label className="h3 m-0">
          {t(`${this.props.parameter.label}.label`)}
        </Label>
        {this.renderInput()}
      </div>
    );
  }
}

const TranslatedParameter = withTranslation()(Parameter);

export default TranslatedParameter;
