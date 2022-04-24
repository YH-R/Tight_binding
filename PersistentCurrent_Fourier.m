tic

% refer to notes for the theory

% physical quantities
t = 3;      % hopping in eV
q = -1;     % in units of e, electron charge is -e
hbar = 1; 
B = 1e-4;	% 0.86T, moderate magnetic field
% B = 0;	% 8600T, strong magnetic field
a = 1;      % 2.76e-10m, lattice constant

beta = 40;	% inverse of temperature times k_b (eV^-1) (T=290K)
mu = 0;     % chemical potential (eV)

Nx = 1024;     % number of sites in x direction (periodic)
Ny = 32;     % number of sites in y direction (open)

theta = B * a^2 * q / hbar;	% dimensionless magnetic field

gauge = 0.3; % radian

% column vector to store average current
I = zeros([Ny 1]);

for l=1:Nx
    H = blockHamiltonian(t, l, Nx, Ny, theta, gauge); % Hamiltonian (block)
    %disp(H)
    
    S = blockSin(l, Nx, Ny, theta, gauge); % sin block, to calculate current
    
    % Diagonalization
    [P,D] = eig(H);
    
    P = abs(P).^2; % element wise absolute square
    
    d = diag(D); % column vector containing eigenvalues
    %disp(d)
    d = fermi(d, beta, mu); % element wise application of fermi function
    %disp(d)
    
    % supposed to sum over expected value of stationary states,
    % weighted by fermi function of the energy
    % S*P is the expected value, d is the fermi weights, summed over l
    I = I + 2*q/hbar*t/Nx * S * P * d;  
end

ny = 1:Ny; % x values (used for plotting only)
I = I * 0.00024341348; % convert to SI units (ampere)

% Plotting
figure()
hold on;

%set(gca,'fontsize',16);
set(gca,'fontname','times');
%set(gca,'linewidth',1.5);
ylabel('Current (A)');
xlabel('{\it n_y}') ;
plot (ny(1:Ny),I(1:Ny).^2,'DisplayName',"current");

hold off;

toc

% returns 1 block of the Hamiltonian (after Fourier transformed)
function block = blockHamiltonian(t, l, Nx, Ny, theta, gauge)
    block = zeros([Ny, Ny]);
    
    for j=1:Ny-1
        block(j,j+1) = 1; % upper diagonal
        block(j+1,j) = 1; % lower diagonal
    end
    
    for j=1:Ny
        block(j,j) = 2 * cos(2*pi/Nx*l + j*theta + gauge); % diagonal
    end
    
    block = -t * block;
end

% returns the diagonal S matrix, refer to notes 
function sblock = blockSin(l, Nx, Ny, theta, gauge)
    sblock = zeros([Ny, Ny]);
    
    for j=1:Ny
        sblock(j,j) = sin(2*pi/Nx*l + j*theta + gauge); % diagonal
    end
end