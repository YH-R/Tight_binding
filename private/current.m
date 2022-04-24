% current operator
function cur = current(startx, starty, endx, endy, Nx, N, H, P, invP, q, hbar)
    s = flatten(startx, starty, Nx);
    e = flatten(endx, endy, Nx);

    I = zeros(N);
    I(s,e) = -1i * q / hbar * H(s,e);
    I(e,s) = 1i * q / hbar * H(e,s);
    
    sparseI = sparse(I);

    intermediate = sparseI * P; % sparse multiplication
    v = real(sum(invP.' .* intermediate, 1)); % sum(A^T .* B) is faster version of transpose(diag(A*B))
    
    % append information of edge position as last four entries
    cur = horzcat(v, startx, starty, endx, endy);
end
