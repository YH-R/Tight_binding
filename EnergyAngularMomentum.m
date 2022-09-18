% This script plots the Energy vs Angular momentum for each stationary state
% as a scatter plot, with B being the varying parameter

tic

% physical quantities
t = 3;          % units of eV, hopping
hbar = 1;
m = 1;          % mass of electron

Nx = 4;         % number of cells in x direction
Ny = 3;         % number of cells in y direction
N = Nx * Ny;    % total number of cells

B_steps = 100;  % How many percent of one period

% initialize position operators
disp('Initializing X and Y');
X = X_square(Nx, Ny, a, 0);
Y = Y_square(Nx, Ny, a, 0);

E = zeros([N * B_steps, 1]);     % vector to store energies
Lz_Ex = zeros([N * B_steps, 1]);    % vector to store Lz expectations

for j=1:B_steps
    MAX_B = 2*pi*B_steps/100;    % 2 pi for one period of square lattice

    r = MAX_B*j/B_steps;    % dimensionless magnetic field

    H = Hamiltonian_square(t, Nx, Ny, r);

    % Velocity matrices
    Vx = 1/1i/hbar * (X * H - H * X);
    Vy = 1/1i/hbar * (Y * H - H * Y);

    % Angular momentum
    Lz = m * (X * Vy - Y * Vx);

    % Diagonalization
    [P,D] = eig(H);

    d = diag(D); % column vector
    Lzdiag = diag(P' * Lz * P);

    E(1+(j-1)*N:j*N, 1) = d;
    Lz_Ex(1+(j-1)*N:j*N, 1) = Lzdiag;
end

% reshape to group the levels together
E = reshape(reshape(E, N, B_steps)', N*B_steps, 1);
Lz_Ex = reshape(reshape(Lz_Ex, N, B_steps)', N*B_steps, 1);

% Plotting
figure()
scatter(Lz_Ex(:), E(:), 5);
set(gca,'fontsize',16);
set(gca,'fontname','times');
xlabel('$L_z/\hbar$', 'Interpreter', 'latex');
ylabel('E (eV)');
titlestr = strcat('(Nx, Ny)=(',num2str(Nx),',',num2str(Ny),')');
title(titlestr)
