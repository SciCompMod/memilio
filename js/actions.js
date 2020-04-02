function Actions(selector, actions) {
    let $container = $(selector);
    let actionMap = new Map();
    let formatTime = d3.timeFormat("%d.%m.%Y"); // this is actually a function... formatTime(d)
    let listeners = [];

    function addIntervalTo(action) {
        let number_of_intervals = action.intervals.length;
        if(number_of_intervals > 3) {
            return;
        }

        let $intervals = action.node.find('.intervals');
        let now = new Date();
        let start = formatTime(new Date(now.getTime() + 7 * 4 * number_of_intervals * 24 * 60 * 60 * 1000));
        let end = formatTime(new Date(now.getTime() + 7 * 4 * (number_of_intervals + 1) * 24 * 60 * 60 * 1000));

        action.counter++;

        let $interval = $($('#template-time-interval').html());

        $interval.attr('id', action.id + '-' + action.counter);

        $interval
            .find('input.start')
            .val(start);

        $interval
            .find('input.end')
            .val(end);

        action.intervals.push({
            id: action.id + '-' + action.counter,
            start,
            end,
            active: false
        });

        $intervals.append($interval);
    }

    function removeIntervalFrom(action, checkbox) {
        let id = $(checkbox).parent().attr('id');
        $(checkbox).parent().remove();
        action.intervals = action.intervals.filter(i => i.id != id);
    }

    function activateIntervalOn(action, checkbox) {
        let id = $(checkbox).parent().attr('id');
        action.intervals.filter(i => i.id === id)[0].active = true;
    }

    function notify() {
        let data = [];
        actionMap.forEach((v, k, m) => {
            if(v.active) {
                data.push(JSON.parse(JSON.stringify(v)));
            }
        });
        console.log(data);
        listeners.forEach(l => {
            l(data);
        });
    }

    $container
        .on('change', '.action > :checkbox', function() {
            let action = actionMap.get($(this).parent().attr('action'));

            if (!this.checked) {
                action.intervals = [];
                action.node.find('.intervals').empty();
                action.active = false;
            } else {
                action.active = true;
                addIntervalTo(action);
            }
            notify();
        });

    $container
        .on('change', '.intervals :checkbox', function(){
            let action = actionMap.get($(this).closest('.action').attr('action'));
            if (!this.checked) {
                removeIntervalFrom(action, this);
                /*if(action.intervals.length < 3) {
                    addIntervalTo(action); 
                }*/
            } else {
                addIntervalTo(action);
                activateIntervalOn(action, this);
            }
            notify();
        });

    let $template = $($('#template-action').html());
    for(let i =0; i < actions.length; i++) {
        let $el = $template.clone();
        $el.attr('action', 'action-' + i)
        $el.find('.name').html(actions[i]);

        actionMap.set('action-' + i, {
            id: 'action-' + i,
            active: false,
            name: actions[i],
            node: $el,
            intervals: [],
            counter: 1
        });

        $container.append($el);
    }

    return {
        onchange: function(listener) {
            listeners.push(listener);
        }
    }
}