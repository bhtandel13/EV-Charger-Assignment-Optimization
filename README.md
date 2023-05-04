# EV-Charger-Assignment-Optimization

Formulated MINLP to optimize the placement of EV chargers in the East Bay based on current traffic flows and EV adoption rates. All data provided is for Richmond, Berkeley, and Oakland. 

**Google distance matrix and geocoding APIs require a unique key
**Optimization solved with Matlab's Gurobi (10.0.1) interface through cvx

Can test with initial inputs for Richmond (call dem_param.m): 
p1 = [37.9556027777778, -122.4162583333333]; L = 2000; row = 3; col = 5; 

call dem_param.m > pass saved variables to sta_param.m > pass demand and station variables to opt_eval to run optimizaton and assign chargers
