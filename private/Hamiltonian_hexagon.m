% Hamiltonian for hexagon lattice with Peierls substitution
% r = B * q * a^2 / hbar is dimensionless B field
% Landau gauge : (Ax , Ay , Az ) = B (−y , 0 , 0)
function H = Hamiltonian_hexagon(t, Nx, Ny, r)
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
                        H(C,Z) = exp(-1i*3*sqrt(3)/4*y*r); % right
                        H(Z,C) = exp(+1i*3*sqrt(3)/4*y*r);
                    elseif index1==x && index2==y+1 % vertical
                        if rem(x+y,2) == 0
                            H(C,Z) = 1; % up
                            H(Z,C) = 1;
                        % else 
                        % vertical edge is "removed"
                        end
                    % else not nearest neighbours
                    end
                end
            end
        end
    end

    H = -t * H; % hopping
end
