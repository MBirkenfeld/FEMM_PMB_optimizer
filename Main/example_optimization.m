%% example file for the optimization

clc, clear all

lower_boundaries = [0.2];
upper_boundaries = [1.5];
starting_point = [1];

% create function handle with 
handle = @(a, b, c, d) make_sim_example(a, b, c, d); % not pretty, but this is how you have to pass a function handle


% create Optimization_FEMM Object
opti = Optimization_FEMM(handle, lower_boundaries, upper_boundaries,...
        starting_point, 'R', 100, 'air_gap', 2, 'num_bearings', 1,...
        'm_rotor', 200, 'mat_Stator', 'N38', 'mat_Rotor', 'N38');

% optimizating
results = opti.optimize;
save('example.mat', 'results')


















