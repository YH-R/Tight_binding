% X Position operator for square lattice
% x0 is offset, to show that physics does not depend on origin
function X = X_square(Nx, Ny, a, x0)
    N = Nx * Ny;
    X = zeros([N, N]);

    for x = 1:Nx
        for y = 1:Ny
            flat = flatten(x,y,Nx); % flattened indices

            % In position basis,
            % only non-zero values are along the diagonal
            X(flat,flat) = x0 + a * x;
        end
    end
end
