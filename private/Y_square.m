% Y Position operator for square lattice
% y0 is offset
function Y = Y_square(Nx, Ny, a, y0)
    N = Nx * Ny;
    Y = zeros([N, N]);

    for x = 1:Nx
        for y = 1:Ny
            flat = flatten(x,y,Nx); % flattened indices

            % In position basis,
            % only non-zero values are along the diagonal
            Y(flat,flat) = y0 + a * y;
        end
    end
end
