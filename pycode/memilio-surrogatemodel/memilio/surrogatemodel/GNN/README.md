## Explanation of files and order of execution in our process 

1. Data generation : Generate Data 
	1b) with dampings
2. Grid Search : grid search for model architecture for 4 model types
3. ARMAConv extended grid search : add more layers and nodes to ARMAConv models  
4. ARMAConv hyperparameters : chose best ARMAConv model and vary order and iteration in grid search. Later, try share weights and dropout.  
5. ARMAConv Test : test best model architetcure and save model weights
