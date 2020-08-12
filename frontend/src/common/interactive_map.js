import * as d3 from 'd3';
import * as topojson from 'topojson-client';
import {Datasets} from '../redux/app';

class InteractiveMap {
  constructor(node) {
    this.listeners = [];
    this.selected = {
      dataset: null,
      id: null,
      label: null,
    };

    // get actual with and height
    this.width = node.clientWidth;
    this.height = node.clientHeight;

    this.showNames = false;
    this.svg = d3.select(node);
    this.bundeslaender = this.svg.append('g');
    this.landkreise = this.svg.append('g');
    this.tooltip = this.svg.append('text').classed('title', true);

    this.svg.attr('width', this.width).attr('height', this.height);

    const self = this;
    const path = d3.geoPath();

    function zoomTo(d) {
      const [[x0, y0], [x1, y1]] = path.bounds(d);
      return self.svg
        .transition()
        .duration(750)
        .call(
          zoom.transform,
          d3.zoomIdentity
            .translate(self.width / 2, self.height / 2)
            .scale(Math.min(8, 0.9 / Math.max((x1 - x0) / self.width, (y1 - y0) / self.height)))
            .translate(-(x0 + x1) / 2, -(y0 + y1) / 2)
        );
    }

    function zoomed() {
      const {transform} = d3.event;
      self.bundeslaender.attr('transform', transform);
      //self.bundeslaender.attr("stroke-width", 1 / transform.k);
      self.landkreise.attr('transform', transform);
      //self.landkreise.attr("stroke-width", 1 / transform.k);
    }

    function hide(selection) {
      return selection
        .transition()
        .style('opacity', 0)
        .duration(400)
        .on('end', function (d, i) {
          d3.select(this).classed('hidden', true);
        });
    }

    function show(selection) {
      return selection.style('opacity', 0).classed('hidden', false).transition().style('opacity', 1).duration(400);
    }

    const zoom = d3.zoom().scaleExtent([1, 8]).on('zoom', zoomed);

    this.svg
      .on('click', function () {
        hide(self.landkreise.selectAll('path:not(.hidden)'));
        self.bundeslaender.selectAll('.off').classed('off', false);
        zoomTo(self.land).on('end', (_) => {
          self.selected = {
            dataset: null,
            id: null,
            label: null,
          };
          self.notify();
        });
      })
      .call(
        d3
          .zoom()
          .extent([
            [0, 0],
            [300, 570],
          ])
          .scaleExtent([0.1, 3])
          .wheelDelta(() => {
            return -d3.event.deltaY * (d3.event.deltaMode === 1 ? 0.05 : d3.event.deltaMode ? 1 : 0.002);
          })
          .on('zoom', function (d, i) {
            self.bundeslaender.attr('transform', d3.event.transform);
            self.landkreise.attr('transform', d3.event.transform);
          })
      );

    d3.json('assets/de_topo.json')
      .catch((error) => {
        throw error;
      })
      .then((topo) => {
        this.land = topojson.feature(topo, topo.objects.vg2500_lan);
        this.kreis = topojson.feature(topo, topo.objects.vg2500_krs);

        this.bundeslaender
          .selectAll('path')
          .data(this.land.features)
          .enter()
          .append('path')
          .attr('d', path)
          .classed('land', true)
          .on('click', function (d, i) {
            d3.event.stopPropagation();

            let visible_path = self.landkreise.selectAll('path:not(hidden)');
            let hidden_path = self.landkreise.selectAll('path.hidden');

            hide(visible_path);
            zoomTo(d).on('end', (_) => {
              self.selected = {
                dataset: Datasets.STATES,
                id: parseInt(d.properties.RS),
                label: d.properties.GEN,
              };
              self.notify();
            });
            show(hidden_path.filter((a, b) => a.properties.RS.startsWith(d.properties.RS)));

            self.bundeslaender
              .selectAll('.land')
              .classed('off', false)
              .filter((x) => x.properties.RS !== d.properties.RS)
              .classed('off', true);
          })
          .on('mouseenter', (d, i) => {
            this.tooltip.text(`${d.properties.GEN}`);
          })
          .on('mousemove', (d, i) => {
            this.tooltip.attr('x', d3.event.layerX + 15).attr('y', d3.event.layerY + 20);
          })
          .on('mouseleave', (d, i) => {
            this.tooltip.text(null);
          })
          .append('title')
          .text((d) => `${d.properties.GEN}`);

        if (this.showNames) {
          this.bundeslaender
            .selectAll('text')
            .data(this.land.features)
            .enter()
            .append('text')
            .attr('x', function (d) {
              return path.centroid(d)[0];
            })
            .attr('y', function (d) {
              return path.centroid(d)[1];
            })
            .attr('text-anchor', 'middle')
            .attr('font-size', '20px')
            .text((d) => `${d.properties.GEN}`);
        }

        this.landkreise
          .selectAll('path')
          .data(this.kreis.features)
          .enter()
          .append('path')
          .attr('d', path)
          .classed('hidden', true)
          .classed('kreis', true)
          .on('click', function (d, i) {
            d3.event.stopPropagation();

            self.selected = {
              dataset: Datasets.COUNTIES,
              id: parseInt(d.properties.RS),
              label: d.properties.GEN,
            };
            self.notify();
          })
          .on('mouseover', function (d, i) {
            self.tooltip.text((x) => `${d.properties.GEN}`);
          })
          .on('mousemove', function (d, i) {
            self.tooltip.attr('x', d3.event.layerX + 15).attr('y', d3.event.layerY + 20);
          })
          .on('mouseout', function (d, i) {
            self.tooltip.text(null);
          })
          .append('title')
          .text((d) => `${d.properties.GEN}`);

        zoomTo(this.land);
      });
  }

  notify() {
    this.listeners.forEach((l) => {
      l(this.selected);
    });
  }

  onSelect(listener) {
    this.listeners.push(listener);
  }
}

export default InteractiveMap;
