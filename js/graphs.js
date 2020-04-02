export {
    Graphs
}

class Graphs {
    constructor(selector) {
        this.formatTime = d3.timeFormat("%d.%m.%Y"); // this is actually a function... formatTime(d)
        this.data_read = [];

        this.damping = [];

        this.logscale = false;

        this.xScale = null;
        this.yScale = null;
        this.xScaleMini = null;
        this.yScaleMini = null;

        this.selectLine = null;
        this.selectLineTooltip = [];
        this.selectLineCircle = [];

        this.svg_graphs = d3.select(selector);
        this.margin_top = 300;
        this.margin_bottom = 50;
        this.margin_vertical = this.margin_top + this.margin_bottom;
        this.margin_horizontal = 200;
        this.width = this.svg_graphs.attr("width") - this.margin_horizontal;
        this.height = this.svg_graphs.attr("height") - this.margin_vertical;
        this.heightMini = this.margin_top / 4;

        this.coeffs_critical = [0.2, 0.1, 0.05];
        this.colors = ["#73B0FF", "orange", "red", "black"];

        this.colors_actions = ["#3393ff", "#8333ff ", "#770031", "#00bd1d", "#007765"];

        this.overlay_mini = [];

        this.begin_interval_i = 0,
            this.end_interval_i = 60,
            this.begin_interval_x = 0,
            this.end_interval_x = 60,
            this.mouse_down_on_mini = false;
        this.max_fact = 1;

        this.g = [];
    }

    visualize(data, actions) {
        let self = this;

        function onMouseMove() {
            // get coordinates
            var x0_pos = d3.mouse(this)[0];
            var bisectDate = d3.bisector(function(d) { return d.day; }).left;
            var i = bisectDate(self.data_read[1], self.xScale.invert(x0_pos))
            var selectedData = self.data_read[1][i - 1]
            x0_pos = self.xScale.invert(x0_pos + 0.5 * self.margin_horizontal);

            // set bar to move
            self.selectLine
                .style("opacity", 0.25)
                .attr("x", self.xScale(x0_pos));

            // set date to move with bar
            self.selectLineTooltip[0]
                .style("opacity", 1)
                .attr("text-anchor", "left")
                .attr("alignment-baseline", "middle")
                .attr("x", (self.xScale(x0_pos) + 10))
                .attr("y", self.margin_top)
                .html(self.formatTime(selectedData.day));

            // set cases numbers and circles to move with bar
            for (var i = 0; i < self.data_read[1][0].cases.length; i++) {

                self.selectLineTooltip[i + 1]
                    .style("opacity", 1)
                    .attr("text-anchor", "left")
                    .attr("alignment-baseline", "middle")
                    .attr("x", (self.xScale(x0_pos) + 10))
                    .attr("y", self.margin_top + self.yScale(selectedData.cases[i]))
                    .html(parseInt(selectedData.cases[i]));

                self.selectLineCircle[i]
                    .style("opacity", 1)
                    .attr("cx", self.xScale(x0_pos))
                    .attr("cy", self.margin_top + self.yScale(selectedData.cases[i]))
                    .html(selectedData.cases[i]);

            }
        }

        function onMouseDownMini() {
            // get coordinates
            var x0_pos = d3.mouse(this)[0];
            var bisectDate = d3.bisector(function(d) { return d.day; }).left;
            var i = bisectDate(self.data_read[0], self.xScaleMini.invert(x0_pos));

            self.begin_interval_i = i - 1;
            self.begin_interval_x = x0_pos;

            self.overlay_mini.attr("x", x0_pos);

            self.mouse_down_on_mini = true;
        }

        function onMouseMoveMini() {

            if (self.mouse_down_on_mini) {

                self.end_interval_x = d3.mouse(this)[0];

                if (self.end_interval_x - self.begin_interval_x > 0) {
                    self.overlay_mini.attr("width", self.end_interval_x - self.begin_interval_x);
                }
            }

        }

        function onMouseReleaseMini() {

            self.end_interval_x = d3.mouse(this)[0];

            if (self.end_interval_x - self.begin_interval_x > 0) {
                self.overlay_mini.attr("width", self.end_interval_x - self.begin_interval_x);

                var bisectDate = d3.bisector(function(d) { return d.day; }).left;
                var i = bisectDate(self.data_read[0], self.xScaleMini.invert(self.end_interval_x));

                self.end_interval_i = i - 1;

                self.subsample_data();

                self.update_graph(false);
            }

            self.mouse_down_on_mini = false;
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

            if (self.data_chosen != d3.select(this).attr("alt")) {

                d3.select(this)
                    .transition(200)
                    .attr('fill', '#cdcdcd');
            }

        }


        function onButtonRescaleMouseClick() {

            d3.select(".button_rescale_active").attr("class", "button_rescale")

            var activeElem = d3.select(this);

            self.data_chosen = d3.select(this).attr("alt");
            if (activeElem.node().tagName == "rectangle") {
                activeElem.attr("class", "button_rescale_active")
            }


            d3.selectAll(".button_rescale").transition(200)
                .attr('fill', '#cdcdcd');

            self.logscale = true;
            self.update_graph();

        }

        this.svg_graphs.selectAll("*").remove();
        if (this.data_read.length > 0) { // how to delete data?
            delete(this.data_read[0]);
            delete(this.data_read[1]);
            this.data_read = [];
        }


        if (this.logscale) {
            d3.select("#button_log").attr("class", "button_rescale_active")
                .attr("fill", "#73B0FF")
        } else {
            d3.select("#button_lin").attr("class", "button_rescale_active")
                .attr("fill", "#73B0FF")
        }

        // var xScale = d3.scaleBand().range([0, width]).padding(0);
        this.xScale = d3.scaleTime().range([0, this.width]);
        if (this.logscale) {
            this.yScale = d3.scaleLog().range([this.height, 0]);
            this.max_fact = 5;
        } else {
            this.yScale = d3.scaleLinear().range([this.height, 0]);
        }

        // var xScaleMini = d3.scaleBand().range([0, width]).padding(0);
        this.xScaleMini = d3.scaleLinear().range([0, this.width]);
        if (this.logscale) {
            this.yScaleMini = d3.scaleLog().range([this.heightMini, 0]);
        } else {
            this.yScaleMini = d3.scaleLinear().range([this.heightMini, 0]); // inverted.....
        }

        this.svg_graphs.append('g')
            .append("text")
            .attr("transform", "rotate(-90)")
            .attr("class", "ylabel")
            .attr("text-anchor", "middle")
            .attr("y", 0.25 * this.margin_horizontal)
            .attr("x", -this.margin_top - 0.5 * this.height)
            .text("Fälle");

        this.svg_graphs.append('g')
            .append("text")
            .attr("class", "xlabel")
            .attr("text-anchor", "middle")
            .attr("x", 0.5 * this.width + 0.5 * this.margin_horizontal)
            .attr("y", this.margin_top + 0.75 * this.margin_bottom + this.height)
            .text("Datum");

        this.g[0] = this.svg_graphs.append("rect")
            .attr("width", this.width)
            .attr("height", this.height)
            .attr('fill', 'none')
            .attr("transform", "translate(" + 0.5 * this.margin_horizontal + "," + this.margin_top + ")")
            .attr("pointer-events", "all")
            // .on("mouseover", onMouseOver) //Add listener for the mouseover event
            .on("mousemove", onMouseMove) //Add listener for the mouseout event 
            // .on("mouseout", onMouseOut); //Add listener for the mouseout event 


        this.g[1] = this.svg_graphs.append("rect")
            .attr("width", this.width)
            .attr("height", this.heightMini)
            .attr('fill', 'none')
            .attr("transform", "translate(" + 0.5 * this.margin_horizontal + "," + (this.margin_top - 1.5 * this.heightMini) + ")")
            .attr("pointer-events", "all")
            .on("mousemove", onMouseMoveMini) //Add listener for the mouseover event
            .on("mousedown", onMouseDownMini) //Add listener for the down event 
            // .on("mouseleave", onMouseLeaveMini) //Add listener for the down event 
            .on("mouseup", onMouseReleaseMini) //Add listener for the release event 
            // .on("mouseout", onMouseOutMini); //Add listener for the mouseout event 


        this.g[0] = this.svg_graphs.append("g")
            .attr("id", "main_graph")
            // .attr("width", width) // ist wohl nicht noetig.. mouse events werden auch so durch das rect oben gefangen und nicht durch das g welches ein dimensionsloser Container ist
            // .attr("height", height)
            // .attr('fill', 'none')
            .attr("transform", "translate(" + 0.5 * this.margin_horizontal + "," + this.margin_top + ")");


        this.g[1] = this.svg_graphs.append("g")
            .attr("id", "mini_graph")
            .attr("transform", "translate(" + 0.5 * this.margin_horizontal + "," + (this.margin_top - 1.5 * this.heightMini) + ")");

        this.g[2] = this.svg_graphs.append("g")
            .attr("transform", "translate(" + 0.5 * this.margin_horizontal + "," + 10 + ")");


        this.overlay_mini = this.g[1].append("rect")
            .attr("id", "gray_overlay")
            .style("opacity", 0.2)
            .on("mousemove", onMouseMoveMini) //Add listener for the mouseover event
            .on("mousedown", onMouseDownMini) //Add listener for the down event 
            .on("mouseup", onMouseReleaseMini) //Add listener for the release event 
            .attr("x", this.begin_interval_x)
            .attr("width", this.end_interval_x - this.begin_interval_x)
            .attr("height", this.heightMini); // linker ausgegrauter Bereich


        this.data_read.push(data);

        if (this.begin_interval_i == 0 && this.end_interval_i == 60) {
            this.filteredData = data.filter((d) => { return d.day < data[this.end_interval_i].day });
            this.data_read.push(this.filteredData);
        } else {
            this.data_read.push(data.slice(this.begin_interval_i, this.end_interval_i));
        }


        this.xScale.domain(d3.extent(this.data_read[1], function(d) { return d.day; }));
        this.yScale.domain([0, this.max_fact * d3.max(this.data_read[1], function(d) { return d.cases[0]; })]);


        this.g[0].append("g")
            .attr("id", "xaxis")
            .attr("transform", "translate(0," + this.height + ")")
            .call(d3.axisBottom(this.xScale).tickFormat(d3.timeFormat("%d.%m.%Y")).ticks(6));
        // .tickFormat(d3.timeFormat("%Y-%m-%d")))
        // .selectAll("text")
        // .style("text-anchor", "end")
        // .attr("dx", "-.8em")
        // .attr("dy", ".15em")
        // .attr("transform", "rotate(-65)");

        this.g[0].append("g")
            .attr("id", "yaxis")
            .call(d3.axisLeft(this.yScale).ticks(4));


        var update = false;

        if (this.begin_interval_i == 0 && this.end_interval_i == 60) {
            for (var i = 0; i < this.data_read[0][0].cases.length; i++) {

                this.g[0].append("path")
                    .datum(this.filteredData)
                    .attr("fill", "none")
                    .attr("stroke", this.colors[i])
                    .attr("stroke-width", 1.5)
                    .attr("d", d3.line()
                        .x(d => this.xScale(d.day))
                        .y(d => this.yScale(d.cases[i]))
                    )

            }
        } else {
            update = true;
        }

        // define bars in graphic
        // g.selectAll(".bar_cat1")
        //     .data(filteredData)
        //     .enter();
        // .append("rect")
        // .attr("class", "bar_cat1")
        // .on("mouseover", onMouseOver) //Add listener for the mouseover event
        // .on("mouseout", onMouseOut) //Add listener for the mouseout event       
        // .attr("x", function(d) { return xScale(d.day); })
        // .attr("y", function(d) { return yScale(d.cases); })
        // .attr("width", xScale.bandwidth())
        // .attr("height", function(d) { return height - yScale(d.cases); });

        // mini graph
        this.xScaleMini.domain(d3.extent(this.data_read[0].map(function(d) { return d.day; })));
        this.yScaleMini.domain([1, this.max_fact * d3.max(this.data_read[0], function(d) { return d.cases[0]; })]);

        this.overlay_mini.attr("width", this.end_interval_x - this.begin_interval_x);

        this.g[1].append("g")
            .attr("id", "yaxis")
            .call(d3.axisLeft(this.yScaleMini).tickFormat(function(d) {
                return d;
            }).ticks(4));
        // .attr("transform", "translate(0," + (-margin + 10) + ")");


        this.g[1].append("g")
            .attr("transform", "translate(0," + this.heightMini + ")")
            .call(d3.axisBottom(this.xScaleMini)
                .ticks(6)
                .tickFormat(d3.timeFormat("%d.%m.%Y")));


        if (!update) {
            // Add the line
            for (var i = 0; i < this.data_read[0][0].cases.length; i++) {
                this.g[1].append("path")
                    .datum(this.data_read[0])
                    .attr("fill", "none")
                    .attr("stroke", this.colors[i])
                    .attr("stroke-width", 1.5)
                    .attr("d", d3.line()
                        .x(d => this.xScaleMini(d.day))
                        .y(d => this.yScaleMini(d.cases[i]))
                    )
            }
        } else {
            this.update_graph();
        }



        // actions x-axis
        this.g[2].append("g")
            .attr("transform", "translate(0," + 2 * this.heightMini + ")")
            .call(d3.axisBottom(this.xScaleMini)
                .ticks(6)
                .tickFormat(d3.timeFormat("%d.%m.%Y")));


        // cursor on main graph
        this.selectLine = this.svg_graphs.append("rect")
            .style("opacity", 0)
            .attr("width", 2)
            .attr("height", (d) => { return this.height; })
            .attr("x", 0.75 * this.margin_horizontal)
            .attr("y", this.margin_top);


        for (var i = 0; i < this.data_read[0][0].cases.length; i++) {
            this.selectLineCircle[i] = this.svg_graphs.append('g')
                .append('circle')
                .style("opacity", 0)
                .attr('r', 4)
                .style("fill", this.colors[i])
                .attr("cx", 0.75 * this.margin_horizontal)
                .attr("cy", this.margin_top);
            // .append('circle')
            // .style("opacity", 0)
            // .attr('r', 4)
            // .attr('color', colors[i])
            // .attr("text-anchor", "left")
            // .attr("alignment-baseline", "middle")
            // .attr("x", 0.75 * margin)
            // .attr("y", height - 0.5 * margin);
        }

        for (var i = 0; i < this.data_read[0][0].cases.length + 1; i++) {
            this.selectLineTooltip[i] = this.svg_graphs.append('g')
                .append('text')
                .style("opacity", 0)
                .attr("text-anchor", "left")
                .attr("alignment-baseline", "middle")
                .attr("x", 0.5 * this.margin_horizontal)
                .attr("y", this.margin_top);
        }

        this.updateActions(actions);
    }

    update_graph() {

        this.g[0].selectAll("path").remove();
        this.g[0].select("#xaxis").remove();
        this.g[0].select("#yaxis").remove();

        if (this.logscale) {
            this.yScale = d3.scaleLog().range([this.height, 0]);
            this.max_fact = 5;
        } else {
            this.yScale = d3.scaleLinear().range([this.height, 0]);
        }

        if (this.logscale) {
            this.yScaleMini = d3.scaleLog().range([this.heightMini, 0]);
        } else {
            this.yScaleMini = d3.scaleLinear().range([this.heightMini, 0]); // inverted.....
        }

        this.xScale.domain(d3.extent(this.data_read[1], function(d) { return d.day; }));
        this.yScale.domain([0, this.max_fact * d3.max(this.data_read[1], function(d) { return d.cases[0]; })]);

        this.g[0].append("g")
            .attr("id", "xaxis")
            .attr("transform", "translate(0," + this.height + ")")
            .call(d3.axisBottom(this.xScale).tickFormat(d3.timeFormat("%d.%m.%Y")).ticks(6));
        // .tickFormat(d3.timeFormat("%Y-%m-%d")))
        // .selectAll("text")
        // .style("text-anchor", "end")
        // .attr("dx", "-.8em")
        // .attr("dy", ".15em")
        // .attr("transform", "rotate(-65)");

        this.g[0].append("g")
            .attr("id", "yaxis")
            .call(d3.axisLeft(this.yScale).ticks(4));

        for (var i = 0; i < this.data_read[0][0].cases.length; i++) {

            this.g[0].append("path")
                .datum(this.data_read[1])
                .attr("fill", "none")
                .attr("stroke", this.colors[i])
                .attr("stroke-width", 1.5)
                .attr("d", d3.line()
                    .x(d => this.xScale(d.day))
                    .y(d => this.yScale(d.cases[i]))
                )

        }

        this.g[1].selectAll("path").remove();
        this.g[1].select("#yaxis").remove();

        this.g[1].append("g")
            .attr("id", "yaxis")
            .call(d3.axisLeft(this.yScaleMini).tickFormat(function(d) {
                return d;
            }).ticks(4));


        this.yScaleMini.domain([1, this.max_fact * d3.max(this.data_read[0], function(d) { return d.cases[0]; })]);

        // Add the line
        for (var i = 0; i < this.data_read[0][0].cases.length; i++) {
            this.g[1].append("path")
                .datum(this.data_read[0])
                .attr("fill", "none")
                .attr("stroke", this.colors[i])
                .attr("stroke-width", 1.5)
                .attr("d", d3.line()
                    .x(d => this.xScaleMini(d.day))
                    .y(d => this.yScaleMini(d.cases[i]))
                )
        }

    }


    subsample_data() {
        delete this.data_read[1];
        this.data_read[1] = this.data_read[0].slice(this.begin_interval_i, this.end_interval_i);
    }

    updateActions(actions) {

        this.g[2].selectAll("*").remove();

        this.g[2].append("g")
            .attr("transform", "translate(0," + 2 * this.heightMini + ")")
            .call(d3.axisBottom(this.xScaleMini)
                .ticks(6)
                .tickFormat(d3.timeFormat("%d.%m.%Y")));


        var dummy_start_x = 0,
            dummy_end_x = 0;

        actions
            .forEach((action, index_i) => {
                action.intervals
                    .forEach((interval, index_j) => {
                        dummy_start_x = this.xScaleMini(new Date(interval.start.substr(6, 4), interval.start.substr(3, 2), interval.start.substr(0, 2)));
                        dummy_end_x = this.xScaleMini(new Date(interval.end.substr(6, 4), interval.end.substr(3, 2), interval.end.substr(0, 2)));
                        console.log(action);
                        if (interval.active == true) {
                            this.g[2].append("rect")
                                .attr("class", action.id)
                                .attr("x", d => { return dummy_start_x; })
                                .attr("y", d => { return 30 + this.heightMini + 15 * index_i; })
                                .attr("fill", this.colors_actions[parseInt(action.id[7])])
                                .attr("width", dummy_end_x - dummy_start_x)
                                .attr("height", d => { return 10; });
                        }
                    });
            });
    }
}