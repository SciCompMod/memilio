
/**
 * Simple euler integration
 * @param {*} df Right hand side
 * @param {*} y0
 * @param {*} t
 */
function euler(df, y0, t) {
  var y = tf.tensor1d(y0);
  var t_data = t.dataSync();

  var dt = tf.scalar(t_data[1] - t_data[0])
  var count = t.shape[0];

  var y_vals = [];
  for (var i = 0; i < count; ++i) {
    y_vals.push(y.dataSync());

    dy = tf.tensor1d(df(t_data[i], y));
    y = tf.add(y, tf.mul(dt, dy));
  }

  return tf.tensor2d(y_vals).transpose();
}


function integration_test() {

  const x = tf.linspace(0., 3.14159265 * 2, 20);
  const y = tf.cos(x);
  const ref = tf.sin(x);

  function my_func(t, /* y */) {
    return [Math.cos(t), -Math.sin(t)];
  }

  var res = euler(my_func, [0., 1.], x)
  const [sine, cosine] = tf.split(res, 2, 0);

  sine.print();
  cosine.print();
 
 //TODO: plot both curves

  // in theory this must be one
  var result = tf.add(sine.pow(2), cosine.pow(2));
  return result;
}

// Tiny TFJS train / predict example.
async function run() {

  var result = integration_test();



  // Use the model to do inference on a data point the model hasn't seen.
  // Should print approximately 39.
  document.getElementById('micro-out-div').innerText = result.dataSync();
}

run();
