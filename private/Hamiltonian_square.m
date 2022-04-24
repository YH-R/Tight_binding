% Hamiltonian for square lattice with Peierls substitution
% r = B * q * a^2 / hbar is dimensionless B field
% Landau gauge : (Ax , Ay , Az ) = B (−y , 0 , 0)
function H = Hamiltonian_square(t, Nx, Ny, r)
    N = Nx * Ny;
    H = zeros([N, N]);
    
    for index1 = 1:Nx % creation index
        for index2 = 1:Ny % creation index
            for x = 1:Nx % annihilation index
                for y = 1:Ny % annihilation index
                    % flattened indices
                    C = flatten(index1,index2,Nx);
                    Z = flatten(x,y,Nx);

                    if index1==x+1 && index2==y % right
                        H(C,Z) = exp(-1i*y*r);
                    elseif index1==x-1 && index2==y % left
                        H(C,Z) = exp(+1i*y*r);
                    elseif index1==x && index2==y+1 % top
                        H(C,Z) = 1;
                    elseif index1==x && index2==y-1 % bottom
                        H(C,Z) = 1;
                    % else not nearest neighbours remains zero
                    end
                end
            end
        end
    end

    H = -t * H; % hopping
end