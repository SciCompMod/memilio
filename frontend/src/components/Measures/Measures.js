import React, { Component } from 'react';
import { connect } from 'react-redux';
import { withTranslation } from 'react-i18next';
import {
  Button,
  Collapse,
  CustomInput,
  Modal,
  ModalBody,
  ModalFooter,
  ModalHeader
} from 'reactstrap';

import DateRangePicker from 'react-daterange-picker';

import Moment from 'moment';
import { extendMoment } from 'moment-range';

import {
  addInterval,
  editInterval,
  activateMeasure,
  deactivateMeasure
} from '../../redux/measures';

import 'react-daterange-picker/dist/css/react-calendar.css';
import './Measures.scss';

const moment = extendMoment(Moment);

const format = (format) => {
  return (time) => moment(time).format(format);
};

const State = Object.freeze({
  NEW: Symbol('new'),
  EDIT: Symbol('edit'),
  CLOSED: Symbol('closed')
});

class DateRangeModal extends Component {
  constructor(props) {
    super(props);
    this.state = {
      range: this.sanitize(this.props.interval)
    };
  }

  sanitize(range) {
    if (!range) {
      return null;
    }

    return moment.range(new Date(range.start), new Date(range.end));
  }

  onSelect(range) {
    this.setState({ range });
  }

  select() {
    if (this.props.onSelect) {
      const { start, end } = this.state.range;
      this.props.onSelect({
        ...this.props.interval,
        start: start.toDate().getTime(),
        end: end.toDate().getTime()
      });
    }
  }

  cancel() {
    if (this.props.onCancel) {
      this.props.onCancel();
    }
  }

  render() {
    const { t } = this.props;
    return (
      <Modal
        isOpen={true}
        modalTransition={{ timeout: 10 }}
        className="date-range"
      >
        <ModalHeader>{t('daterangemodal.header')}</ModalHeader>
        <ModalBody>
          <DateRangePicker
            onSelect={this.onSelect.bind(this)}
            value={this.state.range}
            numberOfCalendars={2}
          />
        </ModalBody>
        <ModalFooter>
          <Button
            color="success"
            disabled={this.state.range === null}
            onClick={this.select.bind(this)}
          >
            {t('select')}
          </Button>
          <Button color="danger" onClick={this.cancel.bind(this)}>
            {t('cancel')}
          </Button>
        </ModalFooter>
      </Modal>
    );
  }
}

const TranslatedDateRangeModal = withTranslation()(DateRangeModal);

class Measure extends Component {
  constructor(props) {
    super(props);
    this.state = {
      ...this.props.data,
      interval: null,
      modal: State.CLOSED
    };
  }

  toggle(event) {
    console.log(event.currentTarget);
    const state = this.state;
    this.setState({
      active: !state.active
    });

    const { id } = this.props.data;
    if (state.active) {
      this.props.deactivateMeasure({ id });
    } else {
      this.props.activateMeasure({ id });
    }
  }

  toggleModal() {
    const state = this.state;
    this.setState({
      modal: !state.modal
    });
  }

  newInterval() {
    this.setState({
      modal: State.NEW,
      interval: null
    });
  }

  editInterval(interval) {
    this.setState({
      modal: State.EDIT,
      interval
    });
  }

  close() {
    this.setState({
      modal: State.CLOSED,
      interval: null
    });
  }

  onSelect(interval) {
    switch (this.state.modal) {
      case State.NEW:
        this.props.addInterval({
          measure: this.props.data.id,
          ...interval
        });
        break;

      case State.EDIT:
        this.props.editInterval({
          measure: this.props.data.id,
          ...interval
        });
        break;

      default:
        break;
    }

    this.close();
  }

  renderModal() {
    return this.state.modal !== State.CLOSED ? (
      <TranslatedDateRangeModal
        interval={this.state.interval}
        onCancel={() => this.close()}
        onSelect={this.onSelect.bind(this)}
      />
    ) : (
      ''
    );
  }

  render() {
    const m = this.props.data;
    const t = this.props.t;
    const dateformat = t('dateformat.date');
    return (
      <div className="measure">
        <CustomInput
          type="switch"
          id={m.label}
          checked={this.state.active}
          onChange={(event) => this.toggle(event)}
          className="d-inline"
        />
        <span className="h3 ml-0">{t(m.label)}</span>
        {this.state.active ? (
          <Button
            onClick={() => this.newInterval()}
            size="sm"
            color="green"
            className="ml-2 py-0"
          >
            {t('new')}
          </Button>
        ) : (
          <span></span>
        )}
        <Collapse isOpen={this.state.active}>
          <div className="intervals">
            {m.intervals.map((interval, i) => {
              return (
                <div className="interval" key={interval.id}>
                  <span className="h4">
                    {format(dateformat)(interval.start)} -{' '}
                    {format(dateformat)(interval.end)}
                  </span>
                  <Button
                    color="warning"
                    onClick={() => this.editInterval(interval)}
                  >
                    <i className="fa fa-edit"></i>
                  </Button>
                </div>
              );
            })}
          </div>
        </Collapse>
        {this.renderModal()}
      </div>
    );
  }
}

const TranslatedMeasure = connect(null, {
  addInterval,
  editInterval,
  activateMeasure,
  deactivateMeasure
})(withTranslation()(Measure));

class Measures extends Component {
  render() {
    const { t } = this.props;
    return (
      <div className="measures">
        <div className="header">{t('measures.title')}</div>
        <div className="content">
          {this.props.measures.map((m, i) => {
            return <TranslatedMeasure key={i} data={m} />;
          })}
        </div>
      </div>
    );
  }
}

const mapState = (state) => {
  return {
    measures: state.measures
  };
};

const MeasuresTranslated = withTranslation()(Measures);
const MeasuresConnected = connect(mapState, {})(MeasuresTranslated);

export { MeasuresConnected as Measures };
