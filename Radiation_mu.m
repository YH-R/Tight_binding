% this program plots the radiation power and torque against chemical 
% potential at specific magnetic field and temperature

tic

% physical quantities
t = 3; % units of eV, hopping
q = -1; % Units of e, electron charge
hbar = 1; 
e_over_kb = 11603;
e = 1.609e-19;
Nx = 20; % number of cells in x direction
Ny = 20; % number of cells in y direction
N = Nx * Ny; % total number of cells

B = 14000 * 1.1576763e-4; % units of 8643T, magnetic field
a = 1; % units of 2.76e-10m, lattice constant

c = 2.997e8 / 419382 ; % speed of light
alpha = 7.2973525693e-3; % fine structure constant
mu_list = -4*t + (0:6000)/6000*8*t; % units of eV, chemical potential
beta = 39; % units of eV^-1, inverse temperature

% Peierls phase (proportional to B field)
%r = rand * 2*pi;
r = B * q * a^2 / hbar;

% initialize Hamiltonian (zero matrix)
disp('Initializing H');
% choose correct lattice by uncommenting
%H = Hamiltonian_square(t, Nx, Ny, r);
%H = Hamiltonian_triangle(t, Nx, Ny, r);
H = Hamiltonian_hexagon(t, Nx, Ny, r);

% stationary states and energy levels
disp('Diagonalizing H');
[P,D] = eig(H);

%disp("P")
%disp(P)

% Matlab diagonalization produces unitary matrix P already
disp('Inverting P');
%invP = inv(P);
invP = ctranspose(P);


% initialize position operators
disp('Initializing X and Y');
X = zeros(N);
Y = zeros(N);

for x = 1:Nx
    for y = 1:Ny
        flat = flatten(x,y,Nx); % flattened indices
        
        % In position basis,        
        % only non-zero values are along the diagonal
        X(flat,flat) = a * x;
        Y(flat,flat) = a * y;
    end
end


% Velocity matrices
disp('Initializing Vx and Vy');
Vx = 1/1i/hbar * (X * H - H * X);
Vy = 1/1i/hbar * (Y * H - H * Y);
    
% Acceleration operators
disp('Initializing Ax and Ay');
Ax = 1/1i/hbar * (Vx * H - H * Vx);
Ay = 1/1i/hbar * (Vy * H - H * Vy);

% change to energy basis
Vx = P' * Vx * P;
Vy = P' * Vy * P;
Ax = P' * Ax * P;
Ay = P' * Ay * P;


% for power
A2 = abs(Ax).^2 + abs(Ay).^2;
A2 = tril(A2); % lower triangular due to step function
    
% for torque
AxVy = real(conj(Ax) .* Vy);
AxVy = tril(AxVy); % lower triangular due to step function


mu_length = length(mu_list);

power_vec = zeros(1, mu_length);
torque_vec = zeros(1, mu_length);

% for matrix multiplication later
E = diag(D); % column vector


for k = 1:mu_length
    mu = mu_list(k);
       
    % temperature terms, refer to notes
    F = fermi(E, beta, mu); % column vector

    power = 4/3 * alpha / hbar / c^2 * (F'*A2*(1-F));
    power_vec(1, k) = power * 2.454e-4; % in W

    torque = -8/3 * alpha * hbar / c^2 * (F'*AxVy*(1-F));
    torque_vec(1,k) = torque * e; % in Nm
end


% Plotting power
figure()
hold on;
scatter(mu_list, power_vec, 5);
scatter(E, zeros(1, length(E)), 5);
set(gca,'fontsize',16);
set(gca,'fontname','times');
%set(gca,'linewidth',1.5);
ylabel('Power (W)');
xlabel('\mu (eV)') ;
axis([-inf inf 0 inf]);
legend(['(Nx, Ny)=(',num2str(Nx),',',num2str(Ny),'), T=', num2str(e_over_kb/beta), 'K, B=', num2str(B/1.1576763e-4), 'T']);

% Plotting torque
figure()
hold on;
scatter(mu_list, torque_vec, 5);
scatter(E, zeros(1, length(E)), 5);
set(gca,'fontsize',16);
set(gca,'fontname','times');
%set(gca,'linewidth',1.5);
ylabel('Torque (Nm)');
xlabel('\mu (eV)') ;
axis([-12 12 -inf inf]);
legend(['(Nx, Ny)=(',num2str(Nx),',',num2str(Ny),'), T=', num2str(e_over_kb/beta), 'K, B=', num2str(B/1.1576763e-4), 'T']);

%hold off; % comment out to plot on same figure

        
toc