from memilio.surrogatemodel.ode_secir_simple.data_generation import generate_data
import os 
import time 
import numpy as np
import tensorflow as tf
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')  # eigentlich unnötig 

### ODE performance 
input_width = 5 
label_width = 25
num_runs = 1 

import time
times = []
for i in range(1000):
  start = time.time()
  data = generate_data(num_runs, path_data, input_width,
                         label_width, save_data = False)

  end = time.time()
  times.append(end - start)
print(np.asarray(times).mean())


# LSTM performance
secirsimple_model = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_150days')
