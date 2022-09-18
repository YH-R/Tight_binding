% this program plots the radiation power and torque against temperature
% at specific magnetic field and chemical potential

tic

% physical quantities
t = 3; % units of eV, hopping
q = -1; % Units of e, electron charge
hbar = 1;
e_over_kb = 11603;
e = 1.60217663e-19;
Nx = 30; % number of cells in x direction
Ny = 20; % number of cells in y direction
N = Nx * Ny; % total number of cells

B = 1 * 1.1576763e-4; % units of 8643T, magnetic field
a = 1; % units of 2.76e-10m, lattice constant

c = 2.99792458e8 / 419382 ; % speed of light
alpha = 7.2973525693e-3; % fine structure constant
mu_list = [-2 -4 -6]; % units of eV, chemical potential
%beta_list = [0.0001]; % units of eV^-1, inverse temperature
temp_list = 10000 + (0:500)*700; % in K
beta_list = e_over_kb./ temp_list;
beta2_list = beta_list.^2;


% Peierls phase (proportional to B field)
r = B * q * a^2 / hbar;

% initialize Hamiltonian (zero matrix)
disp('Initializing H');
H = Hamiltonian_square(t, Nx, Ny, r);

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
X = X_square(Nx, Ny, a, 0);
Y = Y_square(Nx, Ny, a, 0);

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


beta_length = length(beta_list);
power_vec = zeros(1, beta_length);
torque_vec = zeros(1, beta_length);

% for matrix multiplication later
E = diag(D); % column vector

for j = 1:length(mu_list)
    mu = mu_list(j);

    for k = 1:beta_length
        beta = beta_list(k);

        % temperature terms, refer to notes
        F = fermi(E, beta, mu);

        power = 4/3 * alpha / hbar / c^2 * (F'*A2*(1-F));
        power_vec(1, k) = power * 2.454e-4; % in W

        torque = -8/3 * alpha * hbar / c^2 * (F'*AxVy*(1-F));
        torque_vec(1,k) = torque * e; % in Nm
    end


    % Plotting
    hold on;
    figure(1)
    scatter(temp_list, power_vec, 5);

    hold on;
    figure(2)
    scatter(beta2_list, torque_vec, 5);
end

legendcell = strcat('\mu=', strsplit(num2str(mu_list)),'eV');

% setup power plot
figure(1)
set(gca,'fontsize',16);
set(gca,'fontname','times');
%set(gca,'linewidth',1.5);
ylabel('Power (W)');
xlabel('Temperature (K)') ;
%axis([0 1.7 0 1.7]);
titlestr = strcat('(Nx, Ny)=(',num2str(Nx),',',num2str(Ny),'), B=', ...
    num2str(B/1.1576763e-4), 'T');
title(titlestr)
legend(legendcell);
hold off; % comment out to plot on same figure

% setup torque plot
figure(2)
set(gca,'fontsize',16);
set(gca,'fontname','times');
%set(gca,'linewidth',1.5);
ylabel('Torque (Nm)');
xlabel('\beta^2 (eV^{-2})') ;
%axis([0 inf -inf inf]);
titlestr = strcat('(Nx, Ny)=(',num2str(Nx),',',num2str(Ny),'), B=', ...
    num2str(B/1.1576763e-4), 'T');
title(titlestr)
legend(legendcell);
hold off; % comment out to plot on same figure

toc
