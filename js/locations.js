export {
    Locations
}

class Locations {

    constructor(selector) {
        // get actual with and heigt
        
        this.width = $(selector).width();
        this.height = $(selector).height();
        
        $(selector)
            .attr('width', this.width)
            .attr('height', this.height);

        this.showNames = false;
        this.svg = d3.select(selector);
        this.bundeslaender = this.svg.append("g");
        this.landkreise = this.svg.append("g");
        this.tooltip = this.svg.append("text")
                                .attr("font-size", 20)
                                .attr("text-anchor", "start")
                                .attr("stroke-width", 1);
        
        const self = this;
        const path = d3.geoPath();
        
        function zoomTo(d) {
            const [[x0, y0], [x1, y1]] = path.bounds(d);
            return self.svg.transition()
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
            self.bundeslaender.attr("transform", transform);
            //self.bundeslaender.attr("stroke-width", 1 / transform.k);
            self.landkreise.attr("transform", transform);
            //self.landkreise.attr("stroke-width", 1 / transform.k);
        }
        
        function hide(selection) {
            return selection
                .classed("hidden", true)
                .transition()
                    .style("opacity", 0)
                    .duration(750)
                    .on("end", function(d, i) {
                        d3.select(this).attr("visibility", "hidden");
                    });
        }
        
        function show(selection) {
            return selection
                .attr("visibility", "visible")
                .classed("hidden", false)
                .transition()
                .style("opacity", 1)
                    .duration(750);
        }

        const zoom = d3.zoom()
            .scaleExtent([1, 8])
            .on("zoom", zoomed);

        this.svg
            .on('click', function() {
                hide(self.landkreise.selectAll("path:not(.hidden)"));
                //hide(self.landkreise.selectAll("text:not(.hidden)"));
                zoomTo(self.land);
            })
            .call(d3.zoom()
                .extent([[0, 0], [300, 570]])
                .scaleExtent([0.1, 3])
                .wheelDelta(() => {
                    return -d3.event.deltaY * (d3.event.deltaMode === 1 ? 0.05 : d3.event.deltaMode ? 1 : 0.002);
                })
                .on("zoom", function(d, i) {
                    self.bundeslaender.attr("transform", d3.event.transform);
                    self.landkreise.attr("transform", d3.event.transform);
                })
            );


        d3.json("assets/de_topo.json", (error, topo) => {
            if (error) throw error;
            
            this.land = topojson.feature(topo, topo.objects.vg2500_lan);
            this.kreis = topojson.feature(topo, topo.objects.vg2500_krs);

            this.bundeslaender
                .selectAll("path")
                .data(this.land.features)
                .enter()
                .append("path")
                .attr("d", path)
                .attr("fill", "white")
                .attr("stroke", "black")
                .on("click", function(d, i){
                    d3.event.stopPropagation();
                    let visible_path = self.landkreise.selectAll("path:not(hidden)");
                    //let visible_text = self.landkreise.selectAll("text:not(hidden)");
                    let hidden_path = self.landkreise.selectAll("path.hidden");
                   // let hidden_text = self.landkreise.selectAll("text.hidden");

                    hide(visible_path);
                    //hide(visible_text);
                    zoomTo(d);
                    show(hidden_path.filter((a, b) => a.properties.RS.startsWith(d.properties.RS)));
                    //show(hidden_text.filter((a, b) => a.properties.RS.startsWith(d.properties.RS)));
                })
                .on("mouseenter", function(d, i) {
                    d3.select(this).attr("fill", "#dedede");
                    self.tooltip.text(x => `${d.properties.GEN}`);
                })
                .on("mousemove", function(d, i) {
                    self.tooltip
                        .attr("x", (d3.event.layerX + 15))
                        .attr("y", (d3.event.layerY + 20));
                })
                .on("mouseleave", function(d, i){
                    d3.select(this).attr("fill", "white");
                    self.tooltip.text(null);
                })
                .append("title")
                    .text(d => `${d.properties.GEN}`);
            
            if (this.showNames) {
                this.bundeslaender
                    .selectAll("text")
                    .data(this.land.features)
                    .enter()
                    .append("text")
                    .attr("x", function(d) {
                        return path.centroid(d)[0];
                    })
                    .attr("y", function(d) {
                        return path.centroid(d)[1];
                    })
                    .attr("text-anchor", "middle")
                    .attr("font-size", "20px")
                    .text(d => `${d.properties.GEN}`);
            }

            this.landkreise
                .selectAll("path")
                .data(this.kreis.features)
                .enter()
                .append("path")
                .attr("d", path)
                .attr("fill", "white")
                .attr("stroke", "black")
                .attr("opacity", 0)
                .attr("visibility", 'hidden')
                .classed('hidden', true)
                .on('click', function(d, i) {
                    d3.event.stopPropagation();
                })
                .on("mouseover", function(d, i) {
                    d3.select(this).attr("fill", "#dedede");
                    self.tooltip.text(x => `${d.properties.GEN}`);
                })
                .on("mousemove", function(d, i) {
                    self.tooltip
                        .attr("x", (d3.event.layerX + 15))
                        .attr("y", (d3.event.layerY + 20));
                })
                .on("mouseout", function(d, i){
                    d3.select(this).attr("fill", "white");
                    self.tooltip.text(null);
                })
                .append("title")
                    .text(d => `${d.properties.GEN}`);
            /*
            this.landkreise
                .selectAll("text")
                .data(this.kreis.features)
                .enter()
                .append("text")
                .attr("opacity", 0)
                .attr("visibility", 'hidden')
                .classed('hidden', true)
                .attr("x", function(d) {
                    return path.centroid(d)[0];
                })
                .attr("y", function(d) {
                    return path.centroid(d)[1];
                })
                .attr("text-anchor", "middle")
                .attr("font-size", "8px")
                .text(d => `${d.properties.GEN}`);
            */
            zoomTo(this.land);
        });

        /*
        $.getJSON('/assets/de_topo.json')
            .then(topo => {
                //let t = topojson.feature(topo, topo.objects.vg2500_krs)
                //projection.fitSize([500, 570], topo);
                
                context.beginPath();
                path(topojson.mesh(topo));
                context.stroke();

                //console.log(topo);
                /*this.svg_locations.append("g")
                    .selectAll("path")
                    .data(t.features)
                    .enter()
                    .append("path")
                    .attr("d", path)
                    .attr("fill", "white")
                    .attr("stroke", "black");*/
                //.append("title");
                //.text(d => `${d.properties.name}, ${states.get(d.id.slice(0, 2)).name}`)
            

                /*
                let arcs = json.arcs;
                let geometries = json.objects.vg2500_krs.geometries;

                geometries.forEach(o => {
                    console.log(JSON.stringify(o.arcs));
                    //let paths = o.arcs.map(a => arcs[a]);
                    //console.log(paths);
                    /*this.svg_locations
                    .selectAll('path')
                    .data(paths)
                    .enter()
                    .append('path')
                    .style('fill', 'white')
                    .style('stroke', 'black');*
                });*/

                //console.log(Object.keys(json));
                /*json.features = json.features.filter(f => f.geometry.type === 'MultiPolygon');
                
                this.topojson = json;


                let u = this.svg_locations
                    .selectAll('path')
                    .data(this.geojson.features)
                    .enter()
                    .append('path')
                    .style('fill', 'white')
                    .style('stroke', 'black')
                    .attr('d', geoGenerator); */
            //});
            

        this.loc = [{
            name: "K&ouml;ln",
            population: 1061000,
        }, {
            name: "D&uuml;sseldorf",
            population: 612178,
        }, {
            name: "Aachen",
            population: 245885
        }, {
            name: "Bonn",
            population: 318809
        }];

        this.loc_length = [4, 10, 6, 4];
        this.loc_pos = [
            [200, 300],
            [150, 125],
            [75, 425],
            [225, 475]
        ];

        /*
        for (var i = 0; i < this.loc.length; i++) {
            this.loc[i]['id'] = i;
            this.svg_locations.append("circle")
                .attr("id", "name"+i)
                .attr("class", "button_loc")
                .attr("r", 50)
                .attr("cx", this.loc_pos[i][0])
                .attr("cy", this.loc_pos[i][1])
                .attr("fill", "#cdcdcd")
                .attr("alt", i)
                .on("click", onButtonLocMouseClick) //Add listener for the click event
                .on("mouseover", onButtonMouseOver) //Add listener for the mouseover event 
                .on("mouseout", onButtonMouseOut); //Add listener for the mouseout event 

                this.svg_locations.append("text")
                .attr("x", this.loc_pos[i][0] - 4 * this.loc_length[i])
                .attr("y", this.loc_pos[i][1] + 5)
                .attr("alt", i)
                .on("click", onButtonLocMouseClick) //Add listener for the click event
                .on("mouseover", onTextMouseOver)
                .html(this.loc[i].name)
                .attr("font-weight", 600);

        }*/

        this.selected = this.loc[0];
        /*
        this.svg_locations.select('circle#name0')
                .attr("class", "button_loc_active")
                .attr('fill', '#73B0FF');

        function onButtonLocMouseClick() {

            d3.select(".button_loc_active").attr("class", "button_loc")

            var activeElem = d3.select(this);

            self.selected = self.loc[parseFloat(activeElem.attr("alt"))];

            if (activeElem.node().tagName == "circle") {
                activeElem.attr("class", "button_loc_active")
            }


            d3.selectAll(".button_loc").transition(200)
                .attr('fill', '#cdcdcd');

            self.notify();
        }

        function onButtonMouseOver() {
            d3.select(this)
                .transition(400)
                .attr('fill', '#73B0FF');

            d3.select(this).style("cursor", "pointer");
        }

        function onTextMouseOver() {
            d3.select(this).style("cursor", "pointer");
        }


        function onButtonMouseOut() {

            if (self.selected.id != d3.select(this).attr("alt")) {

                d3.select(this)
                    .transition(200)
                    .attr('fill', '#cdcdcd');
            }

        }
        */
    }
    
    notify() {
        this.listeners
            .forEach(l => {
                l(this.selected);
            });
    }

    getPopulation() {
        return this.selected.population;
    }    

    onselect(listener) {
        this.listeners.push(listener);
    }
}