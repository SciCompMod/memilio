//import { Map } from './maps.js'
import { Actions } from './actions.js'
import { Parameters } from './parameters.js'
import { Graphs } from './graphs.js'
import { simulate_seir, makeSeirParam, Damping } from './seir.js'

function main() {
    // Week number according to the ISO-8601 standard
    Date.prototype.getWeek = function () {
        var target = new Date(this.valueOf());
        var dayNr = (this.getDay() + 6) % 7;
        target.setDate(target.getDate() - dayNr + 3);
        var firstThursday = target.valueOf();
        target.setMonth(0, 1);
        if (target.getDay() != 4) {
            target.setMonth(0, 1 + ((4 - target.getDay()) + 7) % 7);
        }
        return 1 + Math.ceil((firstThursday - target) / 604800000);
    }

    $(function () {
        let $container = $('.actions .content');
        //let locations = new Locations("#locations");
        let parameters = new Parameters($('#parameters'));
        let actions = new Actions($container, [{
            label: "75% Home-Office",
            damping: 0.7
        }, {
            label: "Schulschlie&szlig;ungen",
            damping: 0.8,
        }, {
            label: "Kontaktverbot",
            damping: 0.4,
        }]);

        $('svg#graphs')
            .attr('width', $('svg#graphs').width())
            .attr('height', $('svg#graphs').height());

        let graphs = new Graphs("svg#graphs");

        function simulate() {
            let step_size = 0.1;
            let x = parseInt(1 / step_size, 10);
            let days = parseFloat($("#days").val());
            let p = parameters.getParameters();
            let seir_params = makeSeirParam();

            
            seir_params.a = 1.0 / p.incubation;
            seir_params.b = p.contact_rate;
            seir_params.g = 1 / p.infection;
            seir_params.E0 = p.e0;
            seir_params.N = 1000000; //locations.getPopulation();
        
            // TODO: replace by the actual logic
            let action_damping = actions.getActionsDamping(days);
            if(action_damping == null) {
                action_damping = [{day: 0, damping: 1}];
            }
            
            seir_params.dampings = action_damping.map((v, i) => new Damping(v.day, v.damping));
            
            let data = simulate_seir(0, days, step_size, seir_params);
            
            // select only values of the days 
            Object.keys(data)
                .forEach(key => {
                    data[key] = data[key].filter((v, i) => i % x == 0);
                });

            let startDate = d3.timeParse("%d.%m.%Y")("24.02.2020");
            let result = [];
            // restruct data
            for (let i = 0; i < days; i++) {
                result.push({
                    day: new Date(startDate.getTime() + (i * 24 * 60 * 60 * 1000)),
                    cases: [
                        data.S[i],
                        data.E[i],
                        data.I[i],
                        data.R[i]
                    ]
                })
            }


            graphs.visualize(result, actions.getActions());
        }

        /*locations.onselect((arg) => {
            simulate();
        });*/

        actions.onchange((actions) => {
            simulate();
        });

        parameters.onchange(parameters => {
            simulate();
        });


        $('button.simulate')
            .click(function (event) {
                event.preventDefault();
                event.stopPropagation();

                simulate();
            });

        simulate();
    });
}