% Inverse of flatten
% input: (y-1)*Nx + x, Nx
% output: (x,y)
function [x,y] = unflatten(flattened, Nx)
    [Q,R] = quorem(sym(flattened), sym(Nx));
    if R == 0
        x = Nx;
        y = Q;
    else
        x = R;
        y = Q+1;
    end
end