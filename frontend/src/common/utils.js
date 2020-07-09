import {time} from "@amcharts/amcharts4/core";

export const groupBy = (list, key) => {
  return list.reduce(function (groups, item) {
    const val = item[key];
    groups[val] = groups[val] || [];
    groups[val].push(item);
    return groups;
  }, {});
};

export const reduceBy = (list, key) => {
  return;
};

export const merge = (a, b, key) => {
  if (!a && !b) {
    return [];
  }

  if (!a) {
    return b;
  }

  if (!b) {
    return a;
  }

  let uniqueKeys = null;
  if (Array.isArray(key)) {
    if (key.length === 0) {
      return [];
    }

    if (key.length === 1) {
      key = [key[0], key[0]];
    }

    if (key.length > 2) {
      key = key.slice(0, 2);
    }
  } else {
    key = [key, key];
  }

  const aKeys = a.map((e) => e[key[0]]).sort();
  const bKeys = b.map((e) => e[key[1]]).sort();
  uniqueKeys = new Set([...aKeys, ...bKeys]);

  const merged = [];

  for (let k of uniqueKeys) {
    merged.push({
      ...a.find((e) => e[key[0]] === k),
      ...b.find((e) => e[key[1]] === k)
    });
  }

  return merged;
};

export const sumByKey = (list, key) => {
  return list.reduce((acc, val) => acc + val[key], 0);
};

export const renameKey = (list, a, b) => {
  return list.map((e) => {
    e[b] = e[a];
    delete e[a];
    return e;
  });
};

export const calculateDamping = (measures, base_date, days) => {
  var damping = new Array(days).fill(1);

  measures.forEach((measure, index_i) => {
    measure.intervals.forEach((interval) => {
      let start_date = interval.start;
      let end_date = interval.end;

      let start = Math.floor((start_date - base_date) / (1000 * 60 * 60 * 24));
      let end = Math.min(
        days,
        Math.floor((end_date - base_date) / (1000 * 60 * 60 * 24))
      );

      for (var i = start; i < end; i++) {
        if (measures[index_i].damping < damping[i]) {
          damping[i] = measures[index_i].damping;
        }
      }
    });
  });

  // reduce to day where damping changes
  let reduced = [];
  for (let i = 0; i < days; i++) {
    if (reduced.length === 0) {
      reduced.push({
        day: 0,
        damping: damping[i]
      });
      continue;
    }

    if (damping[i] !== reduced[reduced.length - 1].damping) {
      reduced.push({
        day: i,
        damping: damping[i]
      });
    }
  }

  return reduced;
};

/**
 * Rounds the timestamp down to UTC Midnight.
 * @param timestamp {number}
 * @return {number}
 */
export function roundToUTCMidnight(timestamp) {
  return timestamp - (timestamp % (24 * 60 * 60 * 1000));
}