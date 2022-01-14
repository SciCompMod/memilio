# Tutorial: #
-------------

**1. Make your dataset**

Open the "multiply.py".
Define the number of samples you want (good results with >200k). [samples]
Define the number of parallel runs you want (more ==> faster ==> more cpu power used). [n_runs]
Define the number of dampings which should be used for simulations. [n_dampings]
Save and execute it.
(n_runs=40,samples=300000,n_dampings=3 ==> Calculation time on sc-030083l = ca. 2-10 min)
[for test purposes ca. 300 samples ==> all steps are much faster]

**2. Load and make runs**

Execute gui.py.
"File" => "Add ML-Model" => Select "multiple_dampings_learning.py"
Choose Model: (choose multiple_dampings_learning)
Press "Select Dataset" and select your dataset from 1..
Press "Load dataset"	(Loading time: ca. 5-10 min for 300k samples)
If successfully loaded, "Load dataset" and "Run" are green.
Now you can change the learningrate, epochs, etc..
Afterwards, press "Run".
You can run multiple "runs" with different settings (1 dataset + 1 model => multiple runs).
Just change run_name and other settings you want to evalute (like learningrate, epochs, batch size).
If a run is finished, it will be listed under "Finished runs:".
Double click on the finished run to evaluate the run. 

**3. Import runs**
After a finished run, the results will be saved in the dictionary "./models/".
To use and evalute the run whenever you want, import the run.
"File" => "Import runs" => Select all runs you want.
(Import can take some time)
Afterwards the run is listed under "Finished runs". Double click for evaluation.
