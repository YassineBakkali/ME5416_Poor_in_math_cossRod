clc;
clear all;
close all;

% Get the current working directory
currentFolder = pwd;
relativePath = fullfile(currentFolder, 'cossrod_funcs');
addpath(genpath(relativePath));
relativePath = fullfile(currentFolder, 'math_utils');
addpath(genpath(relativePath));
relativePath = fullfile(currentFolder, 'simulator');
addpath(genpath(relativePath));

L = 1; % Length of the rod [m]
r0 = 0.1; % radius of the rod [m]
start = [0; 0; 0];
nElements = 100; % Number of cross-sections
direction = [1.0; 0.0; 0.0];
normal = [0.0; 1.0; 0.0];

damping_cte = 1.0;

% Material Properties %

E0 = 1e5;               % Young modulus of rubber. [N/m^2]
nu = 0.48;              % Poisson ratio of rubber. [-]
G0 = E0 / (2*(1 + nu)); % Shear modulus of rubber. [N/m^2]

rho_material = 2000;           % Density of rubber        [kg/m^3]  
rho_medium = 1000;      % Density of water          [kg/m^3]

params = init_rod(nElements, direction, normal, L, r0, rho_material, rho_medium, E0, start, G0, [], [], [], []);
dl = L / nElements;

dt = dl*0.05;
final_time = 1;
total_steps = final_time / dt;
coss_rod = cossrod_class(params);
simu = simulator_class(coss_rod, damping_cte, dt);
simu = simu.integrate(final_time,total_steps);

% Assume obj.rod.pos_vects is your 3xn matrix where each column is [x; y; z]
x = simu.rod.pos_vects(1, :);  % X coordinates
y = simu.rod.pos_vects(2, :);  % Y coordinates
z = simu.rod.pos_vects(3, :);  % Z coordinates

% Create a new figure
figure;

% Plot the 3D line
plot3(x, y, z, 'LineWidth', 2);

% Adding labels and title for clarity
xlabel('X Axis');
ylabel('Y Axis');
zlabel('Z Axis');
title('3D Visualization of the Rod Centerline');

% Add grid for better visualization
grid on;

% Optionally, set the aspect ratio to equal to ensure all axes are scaled equally
axis equal;