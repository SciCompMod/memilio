import {euler} from '../common/euler.js';
import * as tf from '@tensorflow/tfjs';

const g = 9.80665;

function freefall_ode(t, y) {
  return [y[1], -g];
}

function freefall_exact(t, v0, y0) {
  return [v0 * t + y0 - 0.5 * g * Math.pow(t, 2), v0 - g * t];
}

test('test euler integration', () => {
  const t = 1; // 1 second simulation
  const y0 = 5; // 5 meters of ground
  const v0 = 0; // no start velocity
  const steps = 10000; // simulation steps

  const x = tf.linspace(0, t, steps);

  const exact = freefall_exact(t, v0, y0);
  const result = euler(freefall_ode, [y0, v0], x);
  let [y, v] = tf.split(result, 2, 0);

  y = y.arraySync()[0].pop();
  v = v.arraySync()[0].pop();

  expect(v).toBeCloseTo(exact[1], 2);
  expect(y).toBeCloseTo(exact[0], 2);
});
