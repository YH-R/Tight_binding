% Fermi-Dirac function
% E should be a (column) vector, the rest scalars
function fermi = fermi(E, beta, mu) 
    % large exponential might cause numerical errors?
    fermi = 1./(exp(beta*(E-mu))+1);
end
