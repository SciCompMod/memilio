import { integration_test, euler } from './euler.js';
import {simulate_seir, makeSeirParam} from './seir.js'

function plot_data(seir_data) {


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
      margin: { t: 0 } 
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

function on_slider_change(value)
{
  var beta = value / 1000.;
  document.getElementById('contact_rate_id').innerHTML = beta.toString();

  update_plot(beta);
}

function main()
{
  var beta_slider = document.getElementById('beta_slider');
  beta_slider.oninput = function() {
    on_slider_change(this.value);
  };

  update_plot(1.5);
}

main();


