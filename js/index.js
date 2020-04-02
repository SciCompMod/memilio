import { Locations } from './locations.js'
import { Actions } from './actions.js'
import { Parameters } from './parameters.js'
import { Graphs } from './graphs.js'
import {simulate_seir, makeSeirParam} from './seir.js'

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
        let locations = new Locations("#locations");
        let parameters = new Parameters($('#parameters'));
        let actions = new Actions($container, [
            "75% Home-Office",
            "Schulschlie&szlig;ungen",
            "Kontaktverbot"
        ]);

        $('svg#graphs')
            .attr('width', $('svg#graphs').width())
            .attr('height', $('svg#graphs').height());

        let graphs = new Graphs("svg#graphs");

        locations.onselect((arg) => {
            graphs.visualize(arg);
        });

        actions.onchange((actions) => {
            graphs.updateActions(actions);
        });

        $('button.simulate')
            .click(function (event) {
                event.preventDefault();
                event.stopPropagation();

                let days = parseFloat($("#days").val());
                let p = parameters.getParameters();
                let seir_params = makeSeirParam();

                seir_params.a = 1 / p.incubation;
                seir_params.b = p.contact_rate;
                seir_params.g = 1 / p.infection;
                seir_params.E0 = p.exposed;
                seir_params.N = locations.getPopulation();

                console.log(simulate_seir(0, days, 0.1, seir_params));
            });
    });
}

main()