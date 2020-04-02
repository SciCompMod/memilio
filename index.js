import Locations from './locations.js'

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

$(function(){
    let $container = $('.actions .content');
    let locations = new Locations("#locations");
    let parameters = new Parameters($('#parameters'));
    let actions = new Actions($container, [
        "HomeOffice", 
        "Schlie&szlig;ung von XYZ", 
        "Schulschlie&szlig;ungen", 
        "Kontaktverbot", 
        "Ausgangssperre"
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

    $('button.refresh')
        .click(function(event) {
            event.preventDefault();
            event.stopPropagation();

            let p = parameters.getParameters();
            console.log(p);
        });
});