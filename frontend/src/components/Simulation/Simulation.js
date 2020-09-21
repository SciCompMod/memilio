import React, {Component} from 'react';
import {withTranslation} from 'react-i18next';
import {connect} from 'react-redux';
import {Label, Input, Button, Form, FormGroup, Col, PopoverBody, UncontrolledPopover} from 'reactstrap';
import DateRangePicker from 'react-daterange-picker';
import {getActiveMeasures} from '../../redux/measures';
import {updateModel} from '../../redux/models';
import {getSelectedData} from '../../redux/app';
import {setData} from '../../redux/seir';

import AgeGroupParameters from '../Parameters/AgeGroupParameters';
import {Parameters} from '../Parameters/Parameters';

import {calculateDamping} from '../../common/utils';
import * as seir from '../../common/seir.js';
import * as moment from 'moment';

import './Simulation.scss';

// load secir library if available
let secir = null;
if (window['secirjs']) {
  window.secirjs().then((f) => (secir = f));
}

/**
 * Convert c++ std::vector into JavaScript array.
 *
 * @param {std::vector} vector
 *
 * @returns {Array} Converted vector
 */
const toArray = (vector) => {
  const arr = [];
  for (var i = 0; i < vector.size(); ++i) {
    arr.push(vector.get(i));
  }
  return arr;
};

/**
 * Create new c++ std::vector<size_t> from JavaScript array.
 *
 * @param {Array} array
 *
 * @returns {std::vector<size_t>} Created vector
 */
function createSizeTVector(array) {
  const vec = new secir.VectorSizeT();
  for (let i = 0; i < array.length; i++) {
    vec.push_back(array[i]);
  }
  return vec;
}

/**
 * Component to encaples parameter settings for an simulation
 */
class Simulation extends Component {
  constructor(props) {
    super(props);

    this.state = {
      names: [],
      selected: null,
      popoverOpen: false,
      startDate: moment('2020-02-24'),
      days: 200,
    };

    this.handleSubmit = this.handleSubmit.bind(this);
    this.updateParameter = this.updateParameter.bind(this);
    this.onSelect = this.onSelect.bind(this);
    this.simulate = this.simulate.bind(this);
  }

  /**
   * Update state when properties changes.
   * This happens when values from the the redux store are updated.
   *
   * @param {Object} props New properties
   * @param {Object} state Previous state
   *
   * @return {Object} Updated state
   */
  static getDerivedStateFromProps(props, state) {
    const names = Object.keys(props.models);
    const selected = state.selected;
    return {
      ...state,
      names,
      selected: selected || names[1],
    };
  }

  /**
   * Simulation method for the SECIR model
   */
  simulate_secir() {
    if (!secir) return;

    const model = this.model();
    if (model.name !== 'secir') {
      console.error('Wrong model parameters provided!');
      // TODO: error handling
    }

    function getValue(name) {
      const parameter = model.parameters.find((p) => p.name === name);
      return parameter ? parameter['0-4'] : 0;
    }

    const startDate = this.state.startDate;
    const times = new secir.StageTimes();

    times.set_incubation(getValue('incubation')); // R_2 ^ (-1) + R_3 ^ (-1)
    times.set_infectious_mild(getValue('infmild')); // 4 - 14(=R4 ^ (-1))
    times.set_serialinterval(getValue('serint')); // 4 - 4.4 // R_2^(-1)+0.5*R_3^(-1)
    times.set_hospitalized_to_home(getValue('hosp2home')); // 7 - 16(=R5 ^ (-1))
    times.set_home_to_hospitalized(getValue('home2hosp')); // 2.5 - 7(=R6 ^ (-1))
    times.set_hospitalized_to_icu(getValue('hosp2icu')); // 1 - 3.5(=R7 ^ (-1))
    times.set_icu_to_home(getValue('icu2home')); // 5 - 16(=R8 ^ (-1))
    times.set_infectious_asymp(getValue('infasy')); // (= R9 ^ (-1)= R_3 ^ (-1) + 0.5 * R_4 ^ (-1))
    times.set_icu_to_death(getValue('icu2death')); // 3.5 - 7(=R5 ^ (-1))

    const probs = new secir.Probabilities();
    probs.set_asymp_per_infectious(getValue('alpha')); // 0.01 - 0.16
    probs.set_risk_from_symptomatic(getValue('beta')); // 0.05 - 0.5
    probs.set_hospitalized_per_infectious(getValue('rho')); // 0.1 - 0.35
    probs.set_icu_per_hospitalized(getValue('theta')); // 0.15 - 0.4
    probs.set_dead_per_icu(getValue('delta')); // 0.15 - 0.77
    probs.set_infection_from_contact(getValue('infprob'));

    const vec_times = new secir.VectorSecirParamsStageTimes();
    const vec_probs = new secir.VectorSecirParamsProbabilities();

    vec_times.push_back(times);
    vec_probs.push_back(probs);

    const populations = new secir.Populations(createSizeTVector([1, secir.SecirCompartments.SecirCount.value]));

    populations.set_total(this.props.selected.population);
    populations.set(createSizeTVector([0, secir.SecirCompartments.E.value]), 14400);
    populations.set(createSizeTVector([0, secir.SecirCompartments.C.value]), this.props.start.Confirmed);
    populations.set(createSizeTVector([0, secir.SecirCompartments.I.value]), 50);
    populations.set(createSizeTVector([0, secir.SecirCompartments.H.value]), 20);
    populations.set(createSizeTVector([0, secir.SecirCompartments.U.value]), 10);
    populations.set(createSizeTVector([0, secir.SecirCompartments.R.value]), this.props.start.Recovered);
    populations.set(createSizeTVector([0, secir.SecirCompartments.D.value]), 0);
    populations.set_difference_from_total(
      createSizeTVector([0, secir.SecirCompartments.S.value]),
      this.props.selected.population
    );

    const contact_freq_mat = new secir.ContactFrequencyMatrix();
    contact_freq_mat.set_cont_freq(0.5, 0, 0);

    // emulate some mitigations
    let action_damping = calculateDamping(this.props.measures, startDate, this.state.days);

    action_damping.forEach((v) => {
      const d = new secir.Damping(v.day, v.damping);
      contact_freq_mat.add_damping(d, 0, 0);
      d.delete();
    });

    const uncertain_contact_matrix = new secir.UncertainContactMatrix(contact_freq_mat);

    // set the params required for the simulation
    const param = new secir.SecirParams(contact_freq_mat);
    param.times = vec_times;
    param.probabilities = vec_probs;
    param.populations = populations;

    contact_freq_mat.delete();
    uncertain_contact_matrix.delete();
    populations.delete();

    const result = secir.simulate(0, this.state.days, 0.1, param);

    param.delete();

    // copy data to plain javascript arrays
    const data = {
      t: toArray(result.t),
      S: toArray(result.nb_sus),
      E: toArray(result.nb_exp),
      nb_car: toArray(result.nb_car),
      I: toArray(result.nb_inf),
      nb_hosp: toArray(result.nb_hosp),
      nb_icu: toArray(result.nb_icu),
      R: toArray(result.nb_rec),
      nb_dead: toArray(result.nb_dead),
    };

    result.delete();

    const result2 = Object.values(
      Object.keys(data)
        .filter((k) => k !== 't')
        .reduce((acc, k) => {
          data[k].forEach((value, index) => {
            let date = new Date(startDate + index * 24 * 60 * 60 * 1000);
            date.setHours(0);
            date.setMilliseconds(0);
            date.setMinutes(0);
            date.setSeconds(0);

            date = date.getTime();
            if (!acc[date]) {
              acc[date] = {
                date,
              };
            }
            acc[date] = Object.assign(acc[date], {
              [k]: parseInt(value),
            });
          });

          return acc;
        }, {})
    );

    this.props.setData(result2);
  }

  /**
   * Simulation method for the SEIR model
   */
  simulate_seir() {
    const model = this.model();
    if (model.name !== 'seir') {
      console.error('Wrong model parameters provided!');
      // TODO: error handling
    }

    const startDate = this.state.startDate;

    let step_size = 0.1;
    let x = parseInt(1 / step_size, 10);
    let days = this.state.days;
    let p = model.parameters.reduce((acc, p) => {
      acc[p.name] = p.value || p.default;
      return acc;
    }, {});
    let seir_params = seir.makeSeirParam();

    seir_params.a = 1.0 / p.incubation;
    seir_params.b = p.contact_rate;
    seir_params.g = 1 / p.infection;
    seir_params.E0 = this.props.selected.population / 2;
    seir_params.I0 = this.props.start.Confirmed; //p.i0;
    seir_params.R0 = this.props.start.Recovered;
    seir_params.N = this.props.selected.population;

    let action_damping = calculateDamping(this.props.measures, startDate, days);

    seir_params.dampings = action_damping.map((v, i) => new seir.Damping(v.day, v.damping));

    let data = seir.simulate_seir(0, days, step_size, seir_params);

    // select only values of the days
    // TODO: this probably should be done differently
    Object.keys(data).forEach((key) => {
      data[key] = data[key].filter((v, i) => i % x === 0);
    });

    const result = Object.values(
      Object.keys(data)
        .filter((k) => k !== 't')
        .reduce((acc, k) => {
          data[k].forEach((value, index) => {
            let date = new Date(startDate + index * 24 * 60 * 60 * 1000);
            date.setHours(0);
            date.setMilliseconds(0);
            date.setMinutes(0);
            date.setSeconds(0);

            date = date.getTime();
            if (!acc[date]) {
              acc[date] = {
                date,
              };
            }
            acc[date] = Object.assign(acc[date], {
              [k]: parseInt(value),
            });
          });

          return acc;
        }, {})
    );

    this.props.setData(result);
  }

  /**
   * Wrapper method deciding which simulation model to run.
   */
  simulate() {
    switch (this.state.selected) {
      case 'secir':
        this.simulate_secir();
        break;
      case 'seir':
        this.simulate_seir();
        break;
      default:
        break;
    }
  }

  /**
   * Getter for the currently selected simulation model
   */
  model() {
    return this.props.models[this.state.selected];
  }

  /**
   * Event handler updating the state when a new simulation start date is selected.
   *
   * @param {Date} startDate New start date for the simulation
   */
  onSelect(startDate) {
    this.setState({startDate});
  }

  /**
   * Dummy handler to prevent default form submit handler to be called.
   *
   * @param {Event} event
   */
  handleSubmit(event) {
    event.preventDefault();
  }

  /**
   * Update parameter with new values.
   *
   * @param {Object} parameter
   */
  updateParameter(parameter) {
    this.props.updateModel({
      ...parameter,
      model: this.state.selected,
    });
  }

  /**
   * Renders different parameters depending on the selected model.
   */
  renderParameters() {
    const model = this.props.models[this.state.selected];
    switch (this.state.selected) {
      case 'secir':
        return (
          <AgeGroupParameters
            ageGroups={this.props.ageGroups}
            parameters={model.parameters}
            onChange={this.updateParameter}
          />
        );
      case 'seir':
        return <Parameters parameters={model.parameters} onChange={this.updateParameter} />;
      default:
        return;
    }
  }

  render() {
    const {t} = this.props;
    return (
      <div className="simulation">
        <div className="header">Simulation</div>
        <div className="content">
          <Form onSubmit={this.handleSubmit}>
            <FormGroup row className="mx-0">
              <Label for="model" sm={4}>
                Start Datum:
              </Label>
              <Col sm={8} className="text-right">
                <span className="mr-2">{this.state.startDate.format(t('dateformat.date'))}</span>
                <Button id="startDate" type="button">
                  <i className="fas fa-calendar-day"></i>
                </Button>
                <UncontrolledPopover placement="bottom" trigger="legacy" target="startDate">
                  <PopoverBody>
                    <DateRangePicker
                      onSelect={this.onSelect}
                      value={this.state.startDate}
                      maximumDate={new Date()}
                      numberOfCalendars={1}
                      selectionType="single"
                    />
                  </PopoverBody>
                </UncontrolledPopover>
              </Col>
            </FormGroup>
            <FormGroup row className="mx-0">
              <Label for="days" sm={4}>
                Tage:
              </Label>
              <Col sm={8}>
                <Input
                  type="number"
                  value={this.state.days}
                  onChange={(event) => {
                    this.setState({days: parseFloat(event.target.value)});
                  }}
                />
              </Col>
            </FormGroup>
            <FormGroup row className="mx-0">
              <Label for="model" sm={4}>
                Model:
              </Label>
              <Col sm={8}>
                <Input
                  name="model"
                  type="select"
                  value={this.state.selected}
                  onChange={(event) => {
                    const selected = event.target.value;

                    this.setState({
                      selected,
                    });
                  }}
                >
                  {this.state.names.map((m) => {
                    return (
                      <option key={m} values={m}>
                        {t(m)}
                      </option>
                    );
                  })}
                </Input>
              </Col>
            </FormGroup>
            <FormGroup row className="mx-0">
              <Col>
                <Button
                  onClick={this.simulate}
                  size="sm"
                  disabled={this.state.selected === null || this.props.selected === null}
                >
                  Simulate
                </Button>
              </Col>
            </FormGroup>
          </Form>
          {this.renderParameters()}
        </div>
      </div>
    );
  }
}

const mapState = (state) => {
  let start = getSelectedData(state);
  if (start) {
    start = start.start;
  }
  return {
    start,
    models: state.models,
    selected: state.app.selected,
    ageGroups: state.app.ageGroups,
    measures: getActiveMeasures(state.measures),
  };
};

const Translated = withTranslation()(Simulation);
export default connect(mapState, {setData, updateModel})(Translated);
