function [loss,store]=Loss(y, t, N, m, r_l, sigma, J, M, R, idx_u, store)
%LOSS  Manopt cost function used in the inner solver (rlbfgs).
% Correspond to ALM objective (equation 19) in the paper
%
% This function computes the objective value at the current manifold point (y, t).
% It is intended to be used as:
%   problem.cost = @(yt,store) Loss(yt.y, yt.t, ..., store);
%
% Inputs:
%   y, t     : current manifold variables
%              - y: (3N) x (r_l*m) complex matrix
%              - t: (3N+1) x 1 complex vector
%   N, m, r_l: problem parameters (system size, hierarchy layers, rank)
%   sigma    : penalty parameter 
%   J, M, R  : problem data 
%   idx_u    : linear indices selecting entries of C used in the dual objective term
%   store    : Manopt cache (storedb), shared between cost and gradient evaluations
%
% Outputs:
%   loss     : scalar objective value
%   store    : updated cache

    % Compute and cache intermediate matrices (C and G) shared by Loss and dLoss.
    store = prepareCG(y, t, N, m, r_l, sigma, J, M, R, store);
    C=store.CG.C;
    G=store.CG.G;

    % Objective:
    loss = -sum(C(idx_u))+sigma/2*norm(G)^2-norm(M)^2/(2*sigma);
end
