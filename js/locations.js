export {
    Locations
}

class Locations {

    constructor(selector) {
        this.listeners = [];
        this.data_chosen;

        let self = this;
        var svg_locations = d3.select(selector);

        var loc = ["D&uuml;sseldorf", "K&ouml;ln", "Aachen", "Bonn"];
        var loc_length = [10, 4, 6, 4]
        var loc_ref = ["data_01", "data_02", "data_03", "data_04"];
        var loc_pos = [
            [150, 125],
            [200, 300],
            [75, 425],
            [225, 475]
        ];


        for (var i = 0; i < loc.length; i++) {
            svg_locations.append("circle")
                .attr("class", "button_loc")
                .attr("r", 50)
                .attr("cx", loc_pos[i][0])
                .attr("cy", loc_pos[i][1])
                .attr("fill", "#cdcdcd")
                .attr("alt", loc_ref[i])
                .on("click", onButtonLocMouseClick) //Add listener for the click event
                .on("mouseover", onButtonMouseOver) //Add listener for the mouseover event 
                .on("mouseout", onButtonMouseOut); //Add listener for the mouseout event 

            svg_locations.append("text")
                .attr("x", loc_pos[i][0] - 4 * loc_length[i])
                .attr("y", loc_pos[i][1] + 5)
                .attr("alt", loc_ref[i])
                .on("click", onButtonLocMouseClick) //Add listener for the click event
                .on("mouseover", onTextMouseOver)
                .html(loc[i])
                .attr("font-weight", 600);
        }
    
        function onButtonLocMouseClick() {

            d3.select(".button_loc_active").attr("class", "button_loc")

            var activeElem = d3.select(this);

            //read_and_visualize(activeElem.attr("alt") + ".csv");
            self.select(activeElem.attr("alt") + ".csv");

            self.data_chosen = d3.select(this).attr("alt");
            if (activeElem.node().tagName == "circle") {
                activeElem.attr("class", "button_loc_active")
            }


            d3.selectAll(".button_loc").transition(200)
                .attr('fill', '#cdcdcd');

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
    }
    
    select(arg) {
        this.listeners
            .forEach(l => l(arg))
    }

    onselect(listener) {
        this.listeners.push(listener);
    }
}