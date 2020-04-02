function Actions($container, actions) {
    function addIntervalTo($action) {
        let number_of_intervals = $action.children().length;
        if(number_of_intervals > 3) {
            return;
        }
        let now = new Date();
        let $interval = $($('#template-time-interval').html());
        $interval
            .find(":checkbox")
            .change(function(){
                if(!this.checked) {
                    $interval.remove();
                } else {
                    addIntervalTo($action);
                }
            });

        $interval
            .find('input.start')
            .val(formatTime(new Date(now.getTime() + 7 * 4 * number_of_intervals * 24 * 60 * 60 * 1000)));

        $interval
            .find('input.end')
            .val(formatTime(new Date(now.getTime() + 7 * 4 * (number_of_intervals + 1) * 24 * 60 * 60 * 1000)))

        $action.append($interval);
    }

    let $template = $($('#template-action').html());

    for(let i =0; i < actions.length; i++) {
        let $el = $template.clone();
        $el.find('.name').html(actions[i]);

        let $intervals = $el.find('.intervals')
        $el.find(':checkbox')
            .change(function(){
                if (!this.checked) {
                    $intervals.empty();
                } else {
                    addIntervalTo($intervals);
                }
            });

        $container.append($el);
    }
}