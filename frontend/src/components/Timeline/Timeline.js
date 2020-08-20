import React, {Component} from 'react';
import {connect} from 'react-redux';
import {setCurrentDate} from '../../redux/time';

import * as am4core from '@amcharts/amcharts4/core';

import './Timeline.scss';
import {withTranslation} from 'react-i18next';

/**
 * This Component enables the selection and playback of time. It gets its' start and end time from the store and writes
 * a current time to the store, depending on the slider position.
 */
class Timeline extends Component {
  /** @type Animation */
  sliderAnimation;

  /** @type PlayButton */
  playButton;

  /** @type Slider */
  slider;

  state = {
    /** @type string */
    value: '',
  };

  constructor(props) {
    super(props);
    this.node = React.createRef();
  }

  componentDidMount() {
    //console.log('timeline did mount', this.props.startDate, this.props.endDate);
    /** @type Container */
    const container = am4core.create('timelineDiv', am4core.Container);
    container.width = am4core.percent(100);
    container.layout = 'vertical';

    const sliderContainer = container.createChild(am4core.Container);
    sliderContainer.width = am4core.percent(100);
    sliderContainer.padding(0, 50, 0, 50);
    sliderContainer.layout = 'horizontal';

    this.playButton = sliderContainer.createChild(am4core.PlayButton);
    this.playButton.valign = 'middle';
    this.playButton.events.on('toggled', (event) => {
      if (event.target.isActive) {
        this.play();
      } else {
        this.stop();
      }
    });

    this.slider = sliderContainer.createChild(am4core.Slider);
    this.slider.valign = 'middle';
    this.slider.animationDuration = 0;
    this.slider.margin(0, 0, 0, 0);
    this.slider.opacity = 1;
    this.slider.startGrip.background.fill = am4core.color('#6794DC');
    this.slider.background.opacity = 0.5;
    this.slider.background.fill = am4core.color('#6794DC');
    this.slider.marginLeft = 30;

    const dateLabel = container.createChild(am4core.Label);
    dateLabel.text = '';
    dateLabel.align = 'center';

    this.slider.events.on('rangechanged', () => {
      this.setTime();
      dateLabel.text = this.state.value;
    });

    this.slider.startGrip.events.on('drag', () => {
      this.stop();
      this.sliderAnimation.setProgress(this.slider.start);
    });

    this.sliderAnimation = this.slider.animate({property: 'start', to: 1}, 50000, am4core.ease.linear).pause();
    this.sliderAnimation.events.on('animationended', () => {
      this.playButton.isActive = false;
    });
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    console.log('timeline did update', this.props.startDate, this.props.endDate);
    if (this.props.startDate !== prevProps.startDate || this.props.endDate !== prevProps.endDate) {
      this.setTime();
    }
  }

  /**
   * @param {number} date
   * @return {string}
   */
  createDateString(date) {
    return new Date(date).toLocaleDateString({
      year: 'numeric',
      month: '2-digit',
      day: '2-digit',
    });
  }

  setTime() {
    //console.log('timeline set time', this.props.startDate, this.props.endDate);
    const date = this.props.startDate + (this.props.endDate - this.props.startDate) * this.slider.start;

    this.setState({
      //start: this.createDateString(this.props.startDate),
      //end: this.createDateString(this.props.endDate),
      value: this.createDateString(date),
    });

    this.props.setCurrentDate(date);
  }

  play() {
    if (this.slider) {
      if (this.slider.start >= 1) {
        this.slider.start = 0;
        this.sliderAnimation.start();
      }
      this.sliderAnimation.resume();
      this.playButton.isActive = true;
    }
  }

  stop() {
    this.sliderAnimation.pause();
    this.playButton.isActive = false;
  }

  render(ctx) {
    const {t} = this.props;
    return (
      <div>
        <div id="timelineHeader" className="header">
          {t('timeline.title')}
        </div>
        <div id="timelineDiv" style={{height: '55px'}} />
      </div>
    );
  }
}

function mapState(state) {
  const start = state.time.startDate;
  let end = state.time.endDate;

  return {
    startDate: start,
    endDate: end,
  };
}

const TimelineTranslated = withTranslation()(Timeline);
const ConnectedTimeline = connect(mapState, {setCurrentDate})(TimelineTranslated);

export {ConnectedTimeline as Timeline};
