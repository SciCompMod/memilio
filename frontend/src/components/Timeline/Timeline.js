import React, {Component} from 'react';
import {connect} from 'react-redux';
import {setCurrentDate} from "../../redux/time";

import * as am4core from "@amcharts/amcharts4/core";

import './Timeline.scss';

class Timeline extends Component {
    /** @type Animation */
    sliderAnimation;

    /** @type PlayButton */
    playButton;

    /** @type Slider */
    slider;

    state = {
        /** @type string */
        value: ""
    }

    constructor(props) {
        super(props);
        this.node = React.createRef();
    }

    componentDidMount() {
        /** @type Container */
        let container = am4core.create("timelineDiv", am4core.Container);
        container.width = am4core.percent(100);
        container.height = am4core.percent(75);

        const sliderContainer = container.createChild(am4core.Container);

        sliderContainer.width = am4core.percent(100);
        sliderContainer.valign = "bottom";
        sliderContainer.padding(0, 50, 25, 50);
        sliderContainer.layout = "horizontal";
        sliderContainer.height = 50;

        this.playButton = sliderContainer.createChild(am4core.PlayButton);
        this.playButton.valign = "middle";
        this.playButton.events.on("toggled", event => {
            if (event.target.isActive) {
                this.play();
            } else {
                this.stop();
            }
        })

        this.slider = sliderContainer.createChild(am4core.Slider);
        this.slider.valign = "middle";
        this.slider.margin(0, 0, 0, 0);
        this.slider.background.opacity = 0.3;
        this.slider.opacity = 0.7;
        this.slider.background.fill = am4core.color("#ffffff");
        this.slider.marginLeft = 30;
        this.slider.height = 15;
        this.slider.events.on("rangechanged", () => {
            this.setTime();
        });

        this.slider.startGrip.events.on("drag", () => {
            this.stop();
            this.sliderAnimation.setProgress(this.slider.start);
        });

        this.sliderAnimation = this.slider.animate({property: "start", to: 1}, 50000, am4core.ease.linear).pause();
        this.sliderAnimation.events.on("animationended", () => {
            this.playButton.isActive = false;
        });
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (this.props.startDate !== prevProps.startDate || this.props.endDate !== prevProps.endDate) {
            this.setTime();
        }
    }

    /**
     * @param {number} date
     * @return {string}
     */
    createDateString(date) {
        return new Date(date).toLocaleDateString("de-DE", {
            year: "numeric",
            month: "2-digit",
            day: "2-digit"
        });
    }

    setTime() {
        const date = this.props.startDate + ((this.props.endDate - this.props.startDate) * this.slider.start);

        this.setState({
            start: this.createDateString(this.props.startDate),
            end: this.createDateString(this.props.endDate),
            value: this.createDateString(date)
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
        return (
            <div style={{width: "100%", height: "100%"}}>
                <h3>Start Date: {this.state.start}, End Date: {this.state.end}, Current Date: {this.state.value}</h3>
                <div id="timelineDiv" style={{width: "100%", height: "100%"}}/>
            </div>
        );
    }
}

function mapState(state) {
    const start = state.time.startDate;
    let end = state.time.endDate;

    return {
        startDate: start,
        endDate: end
    }
}

const ConnectedTimeline = connect(mapState, { setCurrentDate })(Timeline);

export {ConnectedTimeline as Timeline};
