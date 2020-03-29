import { integration_test } from './euler.js';

function plot_data(curve_data) {

  var xdata = curve_data.x;
  var ydata = curve_data.y;

	var placeholder = document.getElementById('chart1');
  Plotly.newPlot( placeholder,
    [{
      x: xdata,
      y: ydata
    }],
    {
      margin: { t: 0 } 
    } 
  ); 

}



// Tiny TFJS train / predict example.
async function run() {

  var result = integration_test();

  plot_data(result);

}

run();
