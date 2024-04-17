
% Get the current working directory
currentFolder = pwd;
relativePath = fullfile(currentFolder, 'cossrod_funcs');
addpath(genpath(relativePath));
relativePath = fullfile(currentFolder, 'math_utils');
addpath(genpath(relativePath));
relativePath = fullfile(currentFolder, 'simulator');
addpath(genpath(relativePath));

L = 0.5; % Length of the rod [m]
r0 = 0.1; % radius of the rod [m]
g = [0; 0; -9.81];
start = [0; 0; 0];
nElements = 50; % Number of cross-sections
direction = [1.0; 0.0; 0.0];
normal = [0.0; 1.0; 0.0];

% Material Properties %

E0 = 1e5;           % Young modulus of rubber. [N/m^2]
nu = 0.48;            % Poisson ratio of rubber. [-]
G0 = E0 / (2*(1 + nu)); % Shear modulus of rubber. [N/m^2]
d = 0.02;             % Diameter of rod.         [m]
A0 = pi * d^2 / 4;     % Area of cross section    [m^2]
I0 = pi*(d^4)/64;      % Second moment of area    [m^4]
J0 = I0 * 2;            % Polar moment of inertia  [m^4]
rho = 2000;           % Density of rubber        [kg/m^3]  


params = init_rod(nElements, direction, normal, L, r0, rho, E0, start, G0, [], [], [], []);

coss_rod = cossrod_class(params);
simu = simulator_class(coss_rod);
simu.integrate(10,100);