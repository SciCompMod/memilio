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
