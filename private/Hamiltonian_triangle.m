% Hamiltonian for triangular lattice with Peierls substitution
% r = B * q * a^2 / hbar is dimensionless B field
% Landau gauge : (Ax , Ay , Az ) = B (0 , x , 0)
function H = Hamiltonian_triangle(t, Nx, Ny, r)
    N = Nx * Ny;
    H = zeros([N, N]);
    
    for index1 = 1:Nx % creation index
        for index2 = 1:Ny % creation index
            for x = 1:Nx % annihilation index
                for y = 1:Ny % annihilation index
                    % flattened indices
                    C = flatten(index1,index2,Nx);
                    Z = flatten(x,y,Nx);

                    if index1==x+1 && index2==y % horizontal
                        H(C,Z) = 1; % right
                        H(Z,C) = 1;
                    elseif index1==x && index2==y+1 % anti-diagonal
                        H(C,Z) = exp(+1i*sqrt(3)/2*(x+0.25)*r); % top-right
                        H(Z,C) = exp(-1i*sqrt(3)/2*(x+0.25)*r);
                    elseif index1==x-1 && index2==y+1 % diagonal
                        H(C,Z) = exp(+1i*sqrt(3)/2*(x-0.25)*r); % top-left
                        H(Z,C) = exp(-1i*sqrt(3)/2*(x-0.25)*r);
                    % else not nearest neighbours
                    end
                end
            end
        end
    end

    H = -t * H; % hopping
end