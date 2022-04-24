tic

% refer to notes for the theory

% physical quantities
t = 3;      % hopping in eV
q = -1;     % in units of e, electron charge is -e
hbar = 1; 

Nx = 40;     % number of sites in x direction (periodic)
Ny = 40;     % number of sites in y direction (open)
N = Nx*Ny;

gauge = 1;

B_steps = 100;

energies = zeros([N B_steps]); % matrix to store energies

for j=1:B_steps
    theta = 2*pi*j/B_steps; % dimensionless magnetic field
    
    for l=1:Nx
        H = blockHamiltonian(t, l, Nx, Ny, theta, gauge); % Hamiltonian (block)
        %disp(H)
    
        % Diagonalization
        [P,D] = eig(H);
    
        d = diag(D); % column vector
    
        energies((l-1)*Nx+1:l*Nx, j) = d;
    end
end

%writematrix(energies,'energies.csv')

% Plotting
figure()
hold on;

%set(gca,'fontsize',16);
set(gca,'fontname','times');
%set(gca,'linewidth',1.5);
ylabel('Energy (eV)');
xlabel('{\itB} steps') ;
plot(transpose(energies));

hold off;

toc

% returns 1 block of the Hamiltonian (after Fourier transformed)
% see notes
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
