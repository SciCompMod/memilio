import * as am4core from '@amcharts/amcharts4/core';
import * as am4charts from '@amcharts/amcharts4/charts';
import * as am4maps from '@amcharts/amcharts4/maps';

export default class InteractiveHeatMap {
  /** @type MapChart */
  #map = null;

  /** @type MapPolygonSeries */
  #stateChoroplethSeries = null;

  /** @type MapPolygonSeries */
  #countySeries = null;

  /** @type HeatLegend */
  #stateHeatLegend = null;

  /** @type HeatLegend */
  #countyHeatLegend = null;

  /** @type am4core.Label */
  #dataSetLabel = null;

  /** @type IHeatRule */
  #stateHeatRule = null;

  /** @type boolean */
  #seriesHit = false;

  /** @type number */
  selectedState = -1;

  /** @type number */
  selectedCounty = -1;

  /** @callback StateSelectedCallback
   *  @param newState {StateMetaData} */
  onStateSelected = (newState) => {};

  /** @callback CountySelectedCallback
   *  @param newCounty {any} */
  onCountySelected = (newCounty) => {};

  constructor(id) {
    this.#map = am4core.create(id, am4maps.MapChart);
    this.#map.projection = new am4maps.projections.Mercator();
    this.#map.hiddenState.transitionDuration = 100;
    this.#map.maxPanOut = 0.0;

    this.#dataSetLabel = this.#map.createChild(am4core.Label);
    this.#dataSetLabel.align = 'right';
    this.#dataSetLabel.valign = 'bottom';
    this.#dataSetLabel.marginBottom = 20;
    this.#dataSetLabel.marginRight = 20;
    this.#dataSetLabel.text = 'Dataset: RKI';

    this.#stateChoroplethSeries = new am4maps.MapPolygonSeries();
    this.#stateChoroplethSeries.geodataSource.url = 'assets/german-states.geojson';
    this.#stateChoroplethSeries.useGeodata = true;
    this.#map.series.push(this.#stateChoroplethSeries);

    this.#stateHeatRule = {
      property: 'fill',
      target: this.#stateChoroplethSeries.mapPolygons.template,
      min: am4core.color('#DDD', 1),
      max: am4core.color('#F00', 1),
      logarithmic: false,
    };
    this.#stateChoroplethSeries.heatRules.push(this.#stateHeatRule);

    this.#stateHeatLegend = this.#map.createChild(am4charts.HeatLegend);
    this.#stateHeatLegend.series = this.#stateChoroplethSeries;
    this.#stateHeatLegend.height = am4core.percent(90);
    this.#stateHeatLegend.align = 'left';
    this.#stateHeatLegend.valign = 'middle';
    this.#stateHeatLegend.orientation = 'vertical';
    this.#stateHeatLegend.marginLeft = 20;

    const statePolygonTemplate = this.#stateChoroplethSeries.mapPolygons.template;
    statePolygonTemplate.tooltipText = '{name}: {value}';
    statePolygonTemplate.nonScalingStroke = true;
    statePolygonTemplate.applyOnClones = true;

    const selectStateEvent = (e) => {
      this.#seriesHit = true;
      const item = e.target.dataItem.dataContext;
      this.onStateSelected(item);
      this.stateSelected(item);
    };

    statePolygonTemplate.events.on('hit', selectStateEvent);

    const zoomOutEvent = () => {
      this.onStateSelected(null);
      this.onCountySelected(null);
      this.selectedCounty = -1;
      this.stateSelected({id: -1});
    };

    this.#map.events.on('hit', () => {
      if (!this.#seriesHit) {
        zoomOutEvent();
      }
      this.#seriesHit = false;
    });

    this.#map.events.on('zoomlevelchanged', () => {
      if (this.#map.zoomLevel === 1) {
        zoomOutEvent();
      }
    });

    this.#countySeries = new am4maps.MapPolygonSeries();
    this.#countySeries.geodataSource.url = 'assets/counties/lk_germany.geojson';
    this.#countySeries.useGeodata = true;
    this.#countySeries.hidden = true;

    const countyPolygonTemplate = this.#countySeries.mapPolygons.template;
    countyPolygonTemplate.tooltipText = '{GEN}: {value}';
    countyPolygonTemplate.nonScalingStroke = true;
    countyPolygonTemplate.applyOnClones = true;

    const countyHeatRule = {
      property: 'fill',
      target: countyPolygonTemplate,
      min: am4core.color('#DDD', 1),
      max: am4core.color('#F00', 1),
      logarithmic: false,
    };

    this.#countySeries.heatRules.push(countyHeatRule);

    this.#countyHeatLegend = this.#map.createChild(am4charts.HeatLegend);
    this.#countyHeatLegend.series = this.#countySeries;
    this.#countyHeatLegend.height = am4core.percent(90);
    this.#countyHeatLegend.align = 'left';
    this.#countyHeatLegend.valign = 'middle';
    this.#countyHeatLegend.orientation = 'vertical';
    this.#countyHeatLegend.marginLeft = 20;
    this.#countyHeatLegend.hidden = true;

    countyPolygonTemplate.events.on('hit', (e) => {
      this.#seriesHit = true;
      const item = e.target.dataItem.dataContext;
      this.selectedCounty = item !== null ? parseInt(item.RS, 10) : -1;
      this.onCountySelected(item);
    });

    this.#map.series.push(this.#countySeries);

    this.#map.invalidate();
  }

  showCountySeries() {
    this.#stateChoroplethSeries.hide();
    this.#stateHeatLegend.hide();
    this.#countySeries.show();
    this.#countyHeatLegend.show();

    this.#map.invalidate();
  }

  showStateSeries() {
    this.#countySeries.hide();
    this.#countyHeatLegend.hide();
    this.#stateChoroplethSeries.show();
    this.#stateHeatLegend.show();

    this.#map.invalidate();
  }

  setDataSetName(name) {
    this.#dataSetLabel.text = 'Dataset: ' + name;
  }

  /** @param values {Map<number, number>} */
  setStateValues(values) {
    if (values) {
      for (let stateDatum of this.#stateChoroplethSeries.data) {
        stateDatum.value = values.get(stateDatum.id) ?? 0;
      }

      this.#stateChoroplethSeries.invalidateRawData();
    }
  }

  /** @param values {Map<string, number>} */
  setCountyValues(values) {
    if (!values) {
      return;
    }

    for (let countyDatum of this.#countySeries.data) {
      countyDatum.value = values.get(parseInt(countyDatum.RS, 10)) ?? 0;
    }

    this.#countySeries.invalidateRawData();
  }

  /** @private
   *  @param id {number} */
  stateSelected(id) {
    if (id > 0) {
      this.selectedState = id;
    } else {
      // No state is selected.
      this.selectedState = -1;
      this.#map.goHome();
    }
  }
}
