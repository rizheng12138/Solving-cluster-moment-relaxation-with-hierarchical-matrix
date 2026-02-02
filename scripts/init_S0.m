% init_S0.m
% -------------------------------------------------------------------------
% Generate and save an initial point (yt) for the main solver.
%
% Most users do NOT need to modify this script unless they want to
% regenerate initialization for a new N / different parameters.
% -------------------------------------------------------------------------

clear; clc;
this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(this_file));  
addpath(genpath(fullfile(repo_root, 'src')));

% Set experiment options. Default here uses N=64 for quick initialization.
% You can override settings by adding name-value pairs, e.g.:
%   opts = set_default_opts('N', 128, 'r_l', 120);
opts = set_default_opts('N', 64);

% Unpack frequently used parameters for readability.
N   = opts.N;
r_l = opts.r_l;
m   = opts.m;

% Load pre-generated data file. 
load(data_filename(opts));

% Build sparse matrices V and W used by the gradient computation.
[V, W] = build_VW(N, m, r_l);

% -------------------------------------------------------------------------
% Special initialization target S0
% -------------------------------------------------------------------------
% Construct a target matrix S0 (as a vector) to fit with H(yt.y, yt.t, ...).
% Add speye(3*N+1) here because S is required to be PSD in the main problem: 
% shifting by +I moves the target toward the interior of the PSD cone, which
% typically yields a more stable and "feasible-looking" initialization and
% improves early-iteration behavior. It also helps conditioning.
S0 = reshape(J, 3*N+1, 3*N+1) + speye(3*N+1); 
S0 = S0(:);

% -------------------------------------------------------------------------
% Set up a Manopt problem to find yt such that H(yt) approximates S0
% -------------------------------------------------------------------------
% Variable yt consists of:
%   yt.y : complex Euclidean variable of size (3N) x (r_l*m)
%   yt.t : complex Euclidean variable of size (3N+1) x 1
elems.y = euclideancomplexfactory(3*N, r_l*m);
elems.t = euclideancomplexfactory(3*N+1, 1);
manifold = productmanifold(elems);

problemS.M = manifold;

% Cost: normalized squared Frobenius norm of the mismatch between S0 and H(yt)
problemS.cost  = @(yt) norm(S0 - H(yt.y, yt.t, N, m, r_l))^2 / norm(S0)^2;

% Euclidean gradient of the above cost
problemS.egrad = @(yt) grad_S0(yt.y, yt.t, N, m, r_l, V, W, S0);

% Manopt solver settings
optionsS.maxiter   = 100;  % max iterations for rlbfgs in this initialization stage
optionsS.verbosity = 1;    % printing level 

% Run Manopt (rlbfgs). Empty initial point [] lets Manopt choose a default initialization.
[yt, ~] = rlbfgs(problemS, [], optionsS);

% -------------------------------------------------------------------------
% Save initialization to disk
% -------------------------------------------------------------------------
save(s0_filename(opts), 'yt', 'r_l', 'V', 'W', "-v7.3");
