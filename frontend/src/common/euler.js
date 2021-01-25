import * as tf from '@tensorflow/tfjs';

/**
 * Simple euler integration
 * @param {*} df Right hand side
 * @param {*} y0
 * @param {*} t
 */
function euler(df, y0, t) {
  var y = y0.slice();
  var t_data = t.dataSync();

  var dt = t_data[1] - t_data[0];
  var count = t.shape[0];
  var n_dims = y0.length;

  var y_vals = [];
  for (var i = 0; i < count; ++i) {
    y_vals.push(y.slice());

    var dy = df(t_data[i], y);
    for (var dim = 0; dim < n_dims; ++dim) {
      y[dim] = y[dim] + dt * dy[dim];
    }
  }

  return tf.tensor2d(y_vals).transpose();
}

export {euler};
