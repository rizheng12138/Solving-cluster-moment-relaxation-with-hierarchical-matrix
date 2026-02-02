function [V, W] = build_VW(N, m, r_l)
%BUILD_VW Construct sparse matrices V, W used in the solver.
%
% Inputs:
%   N   : system size 
%   m   : hierarchy layers 
%   r_l : rank 
%
% Outputs:
%   V, W used in the solver .

% ---- V construction ----
Urow = setdiff(1:(3*N+1)*3*N, 3*N+1:3*N+1:(3*N+1)*3*N);
Ucol = 1:9*N^2;

nnzV = 9*N^2*(2-1/2^(m-1));
row_idx = zeros(1, nnzV);
col_idx = zeros(1, nnzV);

for l = 1:m
    len  = 3*N/2^(l-1);
    base = 9*N^2*(2-1/2^(l-2));
    for k = 1:2^(l-1)
        I_idx = ((k-1)*len+1):(k*len);
        [I_grid, J_grid] = ndgrid(I_idx, I_idx);
        lin_idx = I_grid(:) + (J_grid(:)-1) * 3 * N;

        s = base + 9*N^2/4^(l-1)*(k-1) + 1;
        t = base + 9*N^2/4^(l-1)*k;

        row_idx(s:t) = Urow(lin_idx);
        col_idx(s:t) = Ucol(lin_idx);
    end
    col_idx(base+1:base+9*N^2/2^(l-1)) = col_idx(base+1:base+9*N^2/2^(l-1)) + 9*N^2*(l-1);
end

V = sparse(row_idx, col_idx, 1, (3*N+1)^2, m*9*N^2);

% ---- W construction ----
[A, B] = ndgrid(1:3*N, 1:r_l*m);
Krow = (A(:)-1) * r_l*m + B(:);
Kcol = (B(:)-1) * (3*N) + A(:);
K = sparse(Krow, Kcol, 1, 3*N * r_l*m, 3*N * r_l*m);

row_idx = zeros(1, 3*N*r_l*m);
col_idx = zeros(1, 3*N*r_l*m);

for l = 1:m
    idx = (r_l*(l-1) + 1):(r_l*l);
    offsets = (0:3*N-1).' * r_l * m;
    [OffsetGrid, IdxGrid] = ndgrid(offsets, idx);
    lin_idx = OffsetGrid(:) + IdxGrid(:);

    s = 3*N*r_l*(l-1) + 1;
    t = 3*N*r_l*l;

    col_idx(s:t) = lin_idx;
    row_idx(s:t) = col_idx(s:t) + 3*N*r_l*m*(l-1);
end

L = sparse(row_idx, col_idx, 1, m*3*N*r_l*m, 3*N*r_l*m);
W = (L*K).';
end
