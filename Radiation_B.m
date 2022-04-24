% this program plots the radiation power and torque against B field
% at specific temperature and chemical potential

tic

% physical quantities
q = -1; % Units of e, electron charge
hbar = 1; 
t = 3; % units of eV, hopping
e_over_kb = 11603;
e = 1.609e-19;
Nx = 5; % number of cells in x direction
Ny = 5; % number of cells in y direction
N = Nx * Ny; % total number of cells

%B_list = 1.1576763e-4 * 200 * (-250 + (0:500)); % square lattice, 2 periods
%B_list = 1.1576763e-4 * 400 * (-250 + (0:500)); % triangle lattice, 1.5 periods
B_list = 1.1576763e-4 * 100 * (-250 + (0:500)); % hexagon lattice, two periods
B_length = length(B_list);

a = 1; % units of 2.76e-10m, lattice constant

c = 2.997e8 / 419382 ; % speed of light
alpha = 7.2973525693e-3; % fine structure constant

mu = -8; % units of eV, chemical potential
beta = 21; % units of eV^-1, inverse temperature


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


power_vec = zeros(1, B_length);
torque_vec = zeros(1, B_length);


for j = 1:B_length
    B = B_list(j);
    
    % Peierls phase (proportional to B field)
    %r = rand * 2*pi;
    r = B * q * a^2 / hbar;
    
    % initialize Hamiltonian (zero matrix)
    %H = Hamiltonian_square(t, Nx, Ny, r);
    %H = Hamiltonian_triangle(t, Nx, Ny, r);
    H = Hamiltonian_hexagon(t, Nx, Ny, r);

    % stationary states and energy levels
    [P,D] = eig(H); % diagonalization

    % Matlab diagonalization produces unitary matrix P already
    %invP = inv(P);
    invP = ctranspose(P); % inverting P effectively
    
    
    % Velocity matrices
    % initialize velocity operators
    Vx = 1/1i/hbar * (X * H - H * X);
    Vy = 1/1i/hbar * (Y * H - H * Y);
    
    % Acceleration operators
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

    % for matrix multiplication later
    E = diag(D); % column vector
    F = fermi(E, beta, mu); % fermi factors
    
    power = 4/3 * alpha / hbar / c^2 * (F'*A2*(1-F));
    power_vec(1, j) = power * 2.454e-4; % in W
    
    torque = -8/3 * alpha / c^2 * (F'*AxVy*(1-F));
    torque_vec(1,j) = torque * e; % in Nm
end


B_list = B_list * 8643; % convert to T


% Plotting power
figure()
scatter(B_list, power_vec, 5);
set(gca,'fontsize',16);
set(gca,'fontname','times');
%set(gca,'linewidth',1.5);
ylabel('Power (W)');
xlabel('B (T)') ;
axis([-inf inf 0 inf]);
legend(['(Nx, Ny)=(',num2str(Nx),',',num2str(Ny),'), T=', num2str(e_over_kb/beta), 'K, \mu=', num2str(mu), 'eV']);
hold on;


% Plotting torque
figure()
scatter(B_list, torque_vec, 5);
set(gca,'fontsize',16);
set(gca,'fontname','times');
%set(gca,'linewidth',1.5);
ylabel('Torque (Nm)');
xlabel('B (T)') ;
%axis([0 1.7 0 1.7]);
legend(['(Nx, Ny)=(',num2str(Nx),',',num2str(Ny),'), T=', num2str(e_over_kb/beta), 'K, \mu=', num2str(mu), 'eV']);
hold on;

%hold off; % comment out to plot on same figure

toc
