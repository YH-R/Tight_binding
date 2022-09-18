% This script plots the Energy vs Angular momentum for each stationary state
% as a scatter plot, with B being the varying parameter

tic

% physical quantities
t = 3;        % units of eV, hopping
q = -1;       % Units of e, electron charge
hbar = 1;

Nx = 5;       % number of cells in x direction
Ny = 5;       % number of cells in y direction
N = Nx * Ny;  % total number of cells

B_list = 1.1576763e-4 * 200 * (-250 + (0:500)); % square lattice, 2 periods

a = 1; % units of 2.76e-10m, lattice constant
