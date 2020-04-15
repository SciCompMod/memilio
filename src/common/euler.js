export {
    euler,
    integration_test
};

/**
 * Simple euler integration
 * @param {*} df Right hand side
 * @param {*} y0
 * @param {*} t
 */
function euler(df, y0, t) {
    var y = y0.slice();
    var t_data = t.dataSync();
  
    var dt = (t_data[1] - t_data[0])
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
  
  
  function integration_test(n) {
  
    const x = tf.linspace(0., 3.14159265 * 2, n);
    const y = tf.cos(x);
    const ref = tf.sin(x);
  
    function my_func(t, /* y */) {
      return [Math.cos(t), -Math.sin(t)];
    }
  
    var res = euler(my_func, [0., 1.], x)
    const [sine, cosine] = tf.split(res, 2, 0);
  
    sine.print();
    cosine.print();
   
    return {
      x: x.arraySync(),
      y: sine.arraySync()[0]
    };
  }
  