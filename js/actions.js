export {
    Actions
}

class Actions {
    constructor(selector, actions) {
        this.$container = $(selector);
        this.actionMap = new Map();
        this.formatTime = d3.timeFormat("%d.%m.%Y"); // this is actually a function... formatTime(d)
        this.listeners = [];

        let self = this;
        this.$container
            .on('change', '.action > :checkbox', function () {
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
            .on('change', '.intervals :checkbox', function () {
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