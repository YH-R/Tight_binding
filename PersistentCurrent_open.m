% This program calculates the eigenvalues and eigenvectors of the tight
% binding Hamiltonian, and uses them to calculate the expectation values of
% the current operator in the energy eigenstates. The results are exported 
% to a csv file.

tic
% dimensions of system

t = 3; % hopping in eV
q = -1; % electron charge is -e
hbar = 1; 
Nx = 10; % number of cells in x direction
Ny = 10; % number of cells in y direction
N = Nx * Ny; % total number of cells

B = 1.1576763e-4; % ~1T, magnetic field
a = 1; % 2.76e-10m, lattice constant

% dimensionless B field
%r = rand * 2*pi;
r = B * q * a^2 / hbar;


% temperature calculations in python

% initialize Hamiltonian (zero matrix)
% choose the correct Hamiltonian by uncommenting
%H = Hamiltonian_square(t, Nx, Ny, r);
H = Hamiltonian_triangle(t, Nx, Ny, r);
%H = Hamiltonian_hexagon(t, Nx, Ny, r);

% stationary states and energy levels
disp('Diagonalizing H');
[P,D] = eig(H);

% calculate once, not in the loop
% Matlab diagonalization produces unitary matrix P already
disp('Inverting P');
%invP = inv(P);
invP = ctranspose(P);
    
% first row of output is the energy levels
output = transpose(diag(D));
% padding to make size consistent later
% -1 to avoid confusion with position indices
output = horzcat(output, -1, -1, -1, -1);


disp('Calculating current operators');
% loop over all nearest neighbour edges
for startx = 1:Nx
    for starty=1:Ny
        if startx~=Nx % cannot start from the rightmost column
            % horizontal edge (going right)
            endx = startx+1;
            endy = starty;
            
            cur = current(startx, starty, endx, endy, Nx, N, H, P, invP, q, hbar);
            output = vertcat(output, cur); % how to preallocate?
        end
        if starty~=Ny % cannot start from top row
            % vertical edge (going up)
            endx = startx;
            endy = starty+1;
            
            cur = current(startx, starty, endx, endy, Nx, N, H, P, invP, q, hbar);
            output = vertcat(output, cur);
        end
        if startx~=Nx && starty~=1 % cannot start at bottom row and rightmost column (see drawing)
            % diagonal edge (going down and right)
            endx = startx+1;
            endy = starty-1;
            
            cur = current(startx, starty, endx, endy, Nx, N, H, P, invP, q, hbar);
            output = vertcat(output, cur);
        end
    end
end

output = double(output);
% edit path to local directory
% edit name to reflect lattice geometry
output_str = strcat('D:/Cloud/OneDrive/OneDrive - National University of Singapore/Notes/NUS/Physics/FYP/Matlab/code/current_', string(Nx), 'x', string(Ny), '_', string(B)', '_triangle.csv');
writematrix(output, output_str);

toc