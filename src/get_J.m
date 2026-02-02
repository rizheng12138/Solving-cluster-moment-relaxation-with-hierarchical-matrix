function J = get_J(N,h)
%GET_J  Construct the objective coefficient matrix J for the primal problem (P).
%
% In the primal problem (P), the objective is linear in the PSD variable M:
%   minimize   <J, M>    (equivalently, J(:)' * S(:) in vectorized form)
%
% This function builds that sparse matrix J and returns it in vectorized form
%
% Inputs:
%   N : system size 
%   h : model parameter
%
% Output:
%   J : vectorized objective coefficient (length (3N+1)^2), stored as a sparse vector

    J = spalloc(3*N + 1, 3*N + 1, 4*N);

    J(3*N, 3) = J(3*N, 3) - 1;
    J(3, 3*N) = J(3, 3*N) - 1;

    for i = 1 : 1 : N-1
        J(3*i,   3*i+3) = J(3*i,   3*i+3) - 1;
        J(3*i+3, 3*i)   = J(3*i+3, 3*i)   - 1;
    end

    for i = 1 : 1 : N
        J(3*i - 2, 3*N+1) = -h;
        J(3*N+1, 3*i - 2) = -h;
    end

    J = 0.5 * J(:);
end
