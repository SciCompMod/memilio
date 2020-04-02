export {
    Actions
}

class Actions {
    constructor(selector, actions) {
        this.$container = $(selector);
        this.actionMap = new Map();
        this.formatTime = d3.timeFormat("%d.%m.%Y"); // this is actually a function... formatTime(d)
        this.listeners = [];

        this.damping_factors = [0.7, 0.8, 0.4];

        this.base_date = new Date(2020, 2, 24);

        let self = this;
        this.$container
            .on('change', '.action > :checkbox', function() {
                let action = self.actionMap.get($(this).parent().attr('action'));

                if (!this.checked) {
                    action.intervals = [];
                    action.node.find('.intervals').empty();
                    action.active = false;
                } else {
                    action.active = true;
                    self.addIntervalTo(action);
                }
                self.notify();
            });

        this.$container
            .on('change', '.intervals :checkbox', function() {
                let action = self.actionMap.get($(this).closest('.action').attr('action'));
                if (!this.checked) {
                    self.removeIntervalFrom(action, this);
                    /*if(action.intervals.length < 3) {
                        addIntervalTo(action); 
                    }*/
                } else {
                    self.addIntervalTo(action);
                    self.activateIntervalOn(action, this);
                }
                self.notify();
            });

        let $template = $($('#template-action').html());

        for (let i = 0; i < actions.length; i++) {
            let $el = $template.clone();
            $el.attr('action', 'action-' + i)
            $el.find('.name').html(actions[i]);

            this.actionMap.set('action-' + i, {
                id: 'action-' + i,
                active: false,
                name: actions[i],
                node: $el,
                intervals: [],
                counter: 1
            });

            this.$container.append($el);
        }

    }

    addIntervalTo(action) {
        let number_of_intervals = action.intervals.length;
        if (number_of_intervals > 3) {
            return;
        }

        let $intervals = action.node.find('.intervals');
        let now = new Date();
        let start = this.formatTime(new Date(now.getTime() + 7 * 4 * number_of_intervals * 24 * 60 * 60 * 1000));
        let end = this.formatTime(new Date(now.getTime() + 7 * 4 * (number_of_intervals + 1) * 24 * 60 * 60 * 1000));

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

    removeIntervalFrom(action, checkbox) {
        let id = $(checkbox).parent().attr('id');
        $(checkbox).parent().remove();
        action.intervals = action.intervals.filter(i => i.id != id);
    }

    activateIntervalOn(action, checkbox) {
        let id = $(checkbox).parent().attr('id');
        action.intervals.filter(i => i.id === id)[0].active = true;
    }

    getCurrentActions(days_simu) {
        let actions = this.actionMap.forEach((v, k, m) => {
            if (v.active) {
                data.push(JSON.parse(JSON.stringify(v)));
            }
        });
        
        var dummy_start_i = 0,
            dummy_end_i = 0;

        var start = 0,
            end = 0;

        var damping = new Array(days_simu).fill(0);

        actions
            .forEach((action, index_i) => {
                console.log(action.name);
                action.intervals
                    .forEach((interval) => {

                        dummy_start_i = new Date(interval.start.substr(6, 4), interval.start.substr(3, 2), interval.start.substr(0, 2));
                        dummy_end_i = new Date(interval.end.substr(6, 4), interval.end.substr(3, 2), interval.end.substr(0, 2));

                        start = Math.floor((dummy_start_i - this.base_date) / (1000 * 60 * 60 * 24));
                        end = Math.max(days_simu, Math.floor((dummy_end_i - this.base_date) / (1000 * 60 * 60 * 24)));

                        for (var i = start; i < end; i++) {

                            if (this.damping_factors[index_i] < damping[i]) {

                                damping[i] = this.damping_factors[index_i];

                            }
                        }
                    });
            });

        return damping;
    }

    notify() {
        let data = [];
        this.actionMap.forEach((v, k, m) => {
            if (v.active) {
                data.push(JSON.parse(JSON.stringify(v)));
            }
        });
        console.log(data);
        this.listeners.forEach(l => {
            l(data);
        });
    }

    onchange(listener) {
        this.listeners.push(listener);
    }
}