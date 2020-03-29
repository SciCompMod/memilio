import { integration_test } from './euler.js';

// Tiny TFJS train / predict example.
async function run() {

  var result = integration_test();



  // Use the model to do inference on a data point the model hasn't seen.
  // Should print approximately 39.
  document.getElementById('micro-out-div').innerText = result.dataSync();
}

run();
