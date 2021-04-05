Replication files for
"The Effect of Social Distancing on the Reach of an Epidemic in Social Networks"
Journal of Economic Interaction and Coordination
https://doi.org/10.1007/s11403-021-00322-9
by Gutin, Hirano, Hwang, Neary, and Toda

Confirmed to run on Matlab R2019b

All files written by Alexis Akira Toda

Main files:

- getNetwork.m		constructs adjacency matrix of network
- SIR_network_sim3.m	simulate SIR model on network given parameters

Simulation files:

- get_figure_table.m	reproduce figures and tables used in paper for large simulations
- setParam.m		set model parameters
- sim_essential.m	model with essential workers
- sim_gradual.m		model with gradual relaxation of social distancing
- sim_large.m		main result with 1000 replications
- sim_twocountry.m	model with two countries

Data files:

- sim_BA1000.mat	simulation results with BA network, 1000 replications
- sim_ERG1000.mat	simulation results with ERG network, 1000 replications
- sim_WS1000.mat	simulation results with WS network, 1000 replications
- sim_essential.mat	simulation results with essential workers
- sim_twocountry.mat	simulation results with two countries

Instructions:
1. To generate results for large simulations, run "sim_large.m". Please change "type" on line 8 to your favorite network. A simulation with 1000 replications will take a day or two, so please be patient.
2. If you do not want to run the simulation, the results are saved in the three .mat files containing the string "1000", so loading them would suffice.
3. Running "get_figure_table.m" will generate the figures and tables for large simulations.
