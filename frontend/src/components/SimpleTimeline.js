import React, {Component} from 'react';
import * as am4core from '@amcharts/amcharts4/core';

/**
 * This Component enables the selection and playback of time.
 */
export default class SimpleTimeline extends Component {
  static defaultProps = {
    start: 0,
    end: 0,
    value: -1,
    onChange: () => {},
  };

  /** @type Animation */
  #sliderAnimation;

  /** @type PlayButton */
  #playButton;

  /** @type Slider */
  #slider;

  state = {
    /** @type string */
    value: '',
  };

  componentDidMount() {
    /** @type Container */
    const container = am4core.create('timelineDiv', am4core.Container);
    container.width = am4core.percent(100);
    container.layout = 'vertical';

    const sliderContainer = container.createChild(am4core.Container);
    sliderContainer.width = am4core.percent(100);
    sliderContainer.padding(0, 50, 0, 50);
    sliderContainer.layout = 'horizontal';

    this.#playButton = sliderContainer.createChild(am4core.PlayButton);
    this.#playButton.valign = 'middle';
    this.#playButton.background.fill = '#aaa';
    this.#playButton.events.on('toggled', (event) => {
      if (event.target.isActive) {
        this.play();
      } else {
        this.stop();
      }
    });

    this.#slider = sliderContainer.createChild(am4core.Slider);
    this.#slider.valign = 'middle';
    this.#slider.animationDuration = 0;
    this.#slider.margin(0, 0, 0, 0);
    this.#slider.opacity = 1;
    this.#slider.startGrip.background.fill = am4core.color('#aaa');
    this.#slider.background.opacity = 0.5;
    this.#slider.background.fill = am4core.color('#aaa');
    this.#slider.marginLeft = 30;

    this.#slider.events.on('rangechanged', () => {
      const step = parseInt(Math.round((this.props.end - this.props.start - 1) * this.#slider.start, 10));
      if (step !== this.lastStep) {
        this.lastStep = step;
        this.props.onChange(step);
      }
    });

    this.#slider.startGrip.events.on('drag', () => {
      this.stop();
      this.#sliderAnimation.setProgress(this.#slider.start);
    });
  }

  componentDidUpdate(prevProps, prevState) {
    if (prevProps.start !== this.props.start || prevProps.end !== this.props.end) {
      this.#sliderAnimation = this.#slider
        .animate({property: 'start', to: 1}, (this.props.end - this.props.start - 1) * 500, am4core.ease.linear)
        .pause();
      this.#sliderAnimation.events.on('animationended', () => {
        this.#playButton.isActive = false;
      });
      if (prevProps.value === -1 && this.props.value > 0) {
        this.#slider.start = this.props.value;
        this.#sliderAnimation.setProgress(this.#slider.start);
      }
    }
  }

  /**
   * @private
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

  /** @private */
  play() {
    if (this.#slider) {
      if (this.#slider.start >= 1) {
        this.#slider.start = 0;
        this.#sliderAnimation.start();
      }
      this.#sliderAnimation.resume();
      this.#playButton.isActive = true;
    }
  }

  /** @private */
  stop() {
    this.#sliderAnimation.pause();
    this.#playButton.isActive = false;
  }

  render() {
    return <div id="timelineDiv" style={{height: '55px'}} />;
  }
}
