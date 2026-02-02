function store = prepareCG(y, t, N, m, r_l, sigma, J, M, R, store)
%PREPARECG  Compute and cache shared intermediates (C and G) for Loss/dLoss.
%
% This helper is called by both Loss and dLoss to avoid redundant computation
% when Manopt evaluates cost and gradient multiple times at the same point
% (e.g., during line-search, quasi-Newton steps, etc.).
%
% Inputs:
%   y, t     : current manifold variables
%              - y: (3N) x (r_l*m) complex matrix
%              - t: (3N+1) x 1 complex vector
%   N, m, r_l: problem parameters (system size, hierarchy layers, rank)
%   sigma    : penalty parameter 
%   J, M, R  : problem data
%   store    : Manopt cache (storedb).
%
% Output:
%   store    : updated cache.

    % Only compute C and G once per point; reuse if already present in cache.
    if ~isfield(store, 'CG')
        C=J-H(y,t,N,m,r_l);

        G = R*C - 1/sigma*M;

        store.CG = struct('C', C, 'G', G);
    end
end
