% This program calculates the thermal average of the total angular
% momentum as a function of temperature. Multiple values of chemical
% potential can be specified as a vector (mu_list) to plot simultaneously.

tic
% dimensions of system

t = 3; % hopping in eV
q = -1; % electron charge is -e
m = 1; % mass of electron
hbar = 1;
kb = 11603;
Nx = 25; % number of cells in x direction
Ny = 25; % number of cells in y direction
N = Nx * Ny; % total number of cells

B = 1.1576763e-4 * 10; % ~1T, magnetic field
a = 1; % 2.76e-10m, lattice constant

temp_list = 10 + (0:3000)*0.5; % in K
beta_list = kb ./ temp_list;
mu_list = [5 -5]; % chemical potential in eV


% Peierls phase (proportional to B field)
r = B * q * a^2 / hbar;

% initialize Hamiltonian (zero matrix)
disp('Initializing H');
H = Hamiltonian_square(t, Nx, Ny, r);


% initialize position operators
disp('Initializing X and Y');
X = X_square(Nx, Ny, a, 0);
Y = Y_square(Nx, Ny, a, 0);

% Velocity matrices
Vx = 1/1i/hbar * (X * H - H * X);
Vy = 1/1i/hbar * (Y * H - H * Y);

% Angular momentum
Lz = m * (X * Vy - Y * Vx);


% stationary states and energy levels
disp('Diagonalizing H');
[P,D] = eig(H);

% calculate once, not in the loop
% Matlab diagonalization produces unitary matrix P already
disp('Inverting P');
%invP = inv(P);
invP = ctranspose(P);

% Energy levels
E = diag(D);

% only diagonal elements in energy basis needed
Lzdiag = diag(P' * Lz * P); % column

beta_length = length(beta_list);
Lz_vec = zeros(1, beta_length); % angular momentum


% Plotting
figure()
hold on;

for j = 1:length(mu_list)
    mu = mu_list(j);
    for k = 1:beta_length
        beta = beta_list(k);
        ferm = fermi(E, beta, mu); % column

        Lz_vec(1, k) = real(Lzdiag.' * ferm);
    end

    % Plotting
    scatter(temp_list, Lz_vec, 5);
end

legendcell = strcat('\mu = ', strsplit(num2str(mu_list)),'eV');

set(gca,'fontsize',16);
set(gca,'fontname','times');
%set(gca,'linewidth',1.5);
ylabel('$L_z/\hbar$', 'Interpreter', 'latex');
xlabel('T (K)') ;
%axis([0 1.7 0 1.7]);
titlestr = strcat('(Nx, Ny)=(',num2str(Nx),',',num2str(Ny),'), B=', ...
    num2str(B/1.1576763e-4), 'T');
title(titlestr)
legend(legendcell);

toc
