% flattens two integer indices into one, 
% for indexing of Hamiltonian
function flattened = flatten(x,y,Nx)
    flattened = (y-1)*Nx + x;
end