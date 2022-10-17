%% Example file for controlling the simulation
% The idea behind this was to create a simulation class that is defined once for
% all solenoid arrangements in general and then controlled by 
% different scripts, which represent the specific magnet arrangements.
% magnet arrangements
%
% The .syntax for objects is used: first, an object s of the
% class Simulation_FEMM() is created. Then individual methods
% (functions) can be applied to this object: e.g. s.make_study()
% s.get_b_on_line()

clc; 
clear all

%% define geometric parameters:
b = 10;
h = 10;
h_hb = 2*h/5;
R = 100;
r = 1;
z = [-1:1:10];

%% Create magnet objects
% Magnet(width, height, x0, y0, mat_name, magnetization direction, group)       

starter1 = Magnet(b, h_hb, 0, 0, 'N38', 0, 1); % The materials must already be read in the FEMM library
starter2 = Magnet(b, h_hb, 0, 0, 'N38', 0, 2);

% Creation of Halbach arrays for rotor and stator:
hb1 = starter1.make_halbach(360, 0, 5);
hb2 = starter2.make_halbach(0, 360, 5);
hb_magnets = [hb1; hb2];

%% Create simulation object with environment, material and magnets.  
% can be either 'planar'(flat) or 'axi' (axially symmetric). 
% t is the depth for 'planar', for 'axi' t does not matter at all

s_hb = Simulation_FEMM('axi', R, hb_magnets); 
% s = Simulation_FEMM('planar', t, magnets);

% Note: Axial symmetric simulation it is slower. 
% If you test and there is no difference in results between a planar and an
% axial symmetric simulation, use planar.

%% or create a simpler way to create a simulation object with a helper function:
% in this case the example is of a partial planar Halbach Array

s_hb3 = make_sim_example(R,[b,h], 'N38', 'N38');

%% run study
results_hb = s_hb3.make_study(z, r, 'hidden', 1); 
% z and r allow arrays to be passed in and will calculate all values with a
% parallel for-loop if parallel computing toolbox is installed. Otherwise
% it will use a serial for-loop which is much slower. 

% results = [z, r, f_z, f_r, b1_max, b2_max, b3_max]
% b1 = abs(B) under bearing (vec B = Flux Density vector)
% b2 = abs(B) to the right of bearing
% b3 = abs(B) above bearing


%% Perform optimization
% For an example of optimization see example_optimization.m