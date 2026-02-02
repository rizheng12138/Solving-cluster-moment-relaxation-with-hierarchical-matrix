function [val,store]=dLoss(y, t, N, m, r_l, sigma, J, M, R, u, V, W, store)
% DLOSS  Euclidean gradient (egrad) of the Manopt objective Loss.
%
% This function is designed to be used as problem.egrad in Manopt:
%   problem.egrad = @(yt,store) dLoss(yt.y, yt.t, ..., store);
%
% Inputs:
%   y, t     : current manifold variables 
%              - y has size (3N) x (r_l*m)
%              - t has size (3N+1) x 1
%   N, m, r_l: problem parameters (system size, hierarchy layers, rank)
%   sigma    : penalty parameter 
%   J, M, R, u: data defining the objective and constraints
%   V, W     : precomputed sparse operators (constructed by build_VW)
%   store    : Manopt cache (storedb). Used to reuse expensive intermediates.
%
% Outputs:
%   val      : struct with fields 'y' and 't' representing Euclidean gradients
%   store    : updated cache.

    % Ensure expensive intermediates (C and G) are computed once per point and cached in store.
    store = prepareCG(y, t, N, m, r_l, sigma, J, M, R, store);
    G = store.CG.G;

    A = -sigma * G' * R - u;

    r = reshape(A * V, 3*N, 3*N*m);

    s = reshape(y.' * r, m * 3*N * r_l * m, 1);

    % Gradient w.r.t. y
    dLdy = reshape(W * s, 3*N, r_l*m);

    % Gradient w.r.t. t
    dLdt = reshape(A, 3*N+1, 3*N+1) * t;

    % Return gradient in Manopt's expected struct format for a product manifold.
    val = struct('y', dLdy, 't', dLdt);
end
