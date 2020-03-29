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

function update_plot(n_steps)
{
  var result = integration_test(n_steps);
  plot_data(result);
}

window.on_slider_change = function(value)
{
  update_plot(value);
}

update_plot(10);
