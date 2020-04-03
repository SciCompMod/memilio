export {
    Locations
}

class Locations {

    constructor(selector) {
        this.listeners = [];
        this.selected = null;


        let self = this;
        this.svg_locations = d3.select(selector);

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

        }

        this.selected = this.loc[0];
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