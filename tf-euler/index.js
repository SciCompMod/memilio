import { integration_test, euler } from './euler.js';
import {simulate_seir, makeSeirParam} from './seir.js'

function plot_data(seir_data) {

  var plot_type = use_log_scale() ? "log" : "linear";

	var placeholder = document.getElementById('chart1');
  Plotly.newPlot( placeholder,
    [
      {
        x: seir_data.t,
        y: seir_data.S,
        name: "susceptible"
      },
      {
        x: seir_data.t,
        y: seir_data.E,
        name: "exposed"
      },
      {
        x: seir_data.t,
        y: seir_data.I,
        name: "infections"
      },
      {
        x: seir_data.t,
        y: seir_data.R,
        name: "immune / recovered"
      }
    ],
    {
      margin: { t: 0 },
      xaxis: {
        title: "#Days since t0"
      },
      yaxis: {
        type: plot_type,
        autorange: true,
        title: 'Number',
      }
    } 
  ); 

}

function update_plot(beta)
{
  var p = makeSeirParam();
  p.b = beta;

  var result = simulate_seir(0., 400., 0.1, p);
  plot_data(result);
}


function beta() {
  var beta_slider = document.getElementById('beta_slider');
  return beta_slider.value / 1000.;
}

function use_log_scale()
{
  var cb_log = document.getElementById('cb_log_id');
  return cb_log.checked;
}

function main()
{
  var beta_slider = document.getElementById('beta_slider');
  beta_slider.oninput = function() {
    var b = beta();
    document.getElementById('contact_rate_id').innerHTML = b.toString();
    update_plot(b);
  };

  var cb_log = document.getElementById('cb_log_id');
  cb_log.onclick = function() {
    update_plot(beta());
    console.log(use_log_scale());
  }

  update_plot(1.5);
}

main();


