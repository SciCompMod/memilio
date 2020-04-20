import { parseTime } from './utils.js';

export { Actions };

class Actions {
  constructor(selector, actions) {
    this.$container = $(selector);
    this.actionMap = new Map();
    this.formatTime = d3.timeFormat('%d.%m.%Y'); // this is actually a function... formatTime(d)
    this.listeners = [];

    // set the base date always to 24.02.2020
    // note: month index starts at 0 in javascript
    this.base_date = new Date(2020, 1, 24);

    let self = this;
    this.$container.on('change', '.action > :checkbox', function () {
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

    this.$container.on('change', '.intervals :checkbox', function () {
      let action = self.actionMap.get(
        $(this).closest('.action').attr('action')
      );
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

    /*
    this.$container
        .on('change', 'input.start, input.end', () => {

        });
*/
    let $template = $($('#template-action').html());

    for (let i = 0; i < actions.length; i++) {
      let $el = $template.clone();
      $el.attr('action', 'action-' + i);
      $el.find('.name').html(actions[i].label);

      this.actionMap.set('action-' + i, {
        id: 'action-' + i,
        active: false,
        name: actions[i].label,
        node: $el,
        intervals: [],
        counter: 1,
        damping: actions[i].damping
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
    let start = this.formatTime(
      new Date(
        now.getTime() + 7 * 4 * number_of_intervals * 24 * 60 * 60 * 1000
      )
    );
    let end = this.formatTime(
      new Date(
        now.getTime() + 7 * 4 * (number_of_intervals + 1) * 24 * 60 * 60 * 1000
      )
    );
    let interval = {
      id: action.id + '-' + action.counter,
      start,
      end,
      active: false
    };

    action.counter++;

    let $interval = $($('#template-time-interval').html());

    $interval.attr('id', interval.id);

    $interval.find('input.start').val(start);

    $interval.find('input.end').val(end);

    const self = this;
    $interval.find('input.start, input.end').change(function (event) {
      event.stopPropagation();
      event.preventDefault();
      let value = $(this).val();

      if ($(this).is('.start')) {
        interval.start = value;
      } else {
        interval.end = value;
      }
      self.notify();
    });

    action.intervals.push(interval);

    $intervals.append($interval);
  }

  removeIntervalFrom(action, checkbox) {
    let id = $(checkbox).parent().attr('id');
    $(checkbox).parent().remove();
    action.intervals = action.intervals.filter((i) => i.id != id);
  }

  activateIntervalOn(action, checkbox) {
    let id = $(checkbox).parent().attr('id');
    action.intervals.filter((i) => i.id === id)[0].active = true;
  }

  getActions() {
    let actions = [];
    this.actionMap.forEach((v, k, m) => {
      if (v.active && v.intervals.some((i) => i.active)) {
        let copy = JSON.parse(JSON.stringify(v));
        copy.intervals = copy.intervals.filter((i) => i.active);
        actions.push(copy);
      }
    });
    return actions.length == 0 ? null : actions;
  }

  getActionsDamping(days) {
    let actions = this.getActions();

    if (actions == null) {
      return null;
    }

    var damping = new Array(days).fill(1);

    actions.forEach((action, index_i) => {
      action.intervals.forEach((interval) => {
        let start_date = parseTime(interval.start);
        let end_date = parseTime(interval.end);

        let start = Math.floor(
          (start_date - this.base_date) / (1000 * 60 * 60 * 24)
        );
        let end = Math.min(
          days,
          Math.floor((end_date - this.base_date) / (1000 * 60 * 60 * 24))
        );

        for (var i = start; i < end; i++) {
          if (actions[index_i].damping < damping[i]) {
            damping[i] = actions[index_i].damping;
          }
        }
      });
    });

    // reduce to day where damping changes
    let reduced = [];
    for (let i = 0; i < days; i++) {
      if (reduced.length == 0) {
        reduced.push({
          day: 0,
          damping: damping[i]
        });
        continue;
      }

      if (damping[i] != reduced[reduced.length - 1].damping) {
        reduced.push({
          day: i,
          damping: damping[i]
        });
      }
    }

    return reduced;
  }

  notify() {
    let data = this.getActions();
    this.listeners.forEach((l) => {
      l(data);
    });
  }

  onchange(listener) {
    this.listeners.push(listener);
  }
}
