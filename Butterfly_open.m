% This code generates the Hofstadter's butterfly spectrum for square,
% triangle and hexagon lattice (but only 1 at a time)
% refer to notes for the theory

tic

% physical quantities
t = 3;      % hopping in eV
q = -1;     % in units of e, electron charge is -e
hbar = 1;

Nx = 7;      % number of sites in x direction (periodic)
Ny = 10;     % number of sites in y direction (open)
N = Nx*Ny;   % dimension of Hamiltonian

B_steps = 100; % How many points to sample in one period

E = zeros([N B_steps]); % matrix to store energies

for j=1:B_steps
    % To define one period, keep only one uncommented
    MAX_B = 2*pi; % square
    %MAX_B = 8*pi/sqrt(3); % triangle
    %MAX_B = 4*pi/(3*sqrt(3)); % hexagon

    r = MAX_B*j/B_steps; % dimensionless magnetic field

    % Hamiltonian
    % uncommment to choose lattice geometry
    H = Hamiltonian_square(t, Nx, Ny, r);
    %H = Hamiltonian_triangle(t, Nx, Ny, r);
    %H = Hamiltonian_hexagon(t, Nx, Ny, r);
    %disp(H)

    % Diagonalization
    [P,D] = eig(H);

    d = diag(D); % column vector

    E(1:N, j) = d;
end

%writematrix(E,'energies.csv')

% Plotting
figure()
hold on;

%set(gca,'fontsize',16);
set(gca,'fontname','times');
%set(gca,'linewidth',1.5);
ylabel('Energy (eV)');
xlabel('{\itB} steps') ;
plot(transpose(E));

hold off;

toc
