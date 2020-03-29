import { integration_test } from './euler.js';

// Tiny TFJS train / predict example.
async function run() {

  var result = integration_test();

  // Now print the resulting array of ones on screen by overwriting the div
  document.getElementById('micro-out-div').innerText = result.dataSync();
}

run();
