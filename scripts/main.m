% main.m
% -------------------------------------------------------------------------
% Entry point for reproducing the benchmark in:
%   "Solving Cluster Moment Relaxation with Hierarchical Matrix"
%
% Typical usage:
%   - Run this file directly to reproduce a default experiment (e.g., N=64).
%   - To change the system size or algorithm parameters, edit the opts line below.
% -------------------------------------------------------------------------
clear; clc; 
this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(this_file));
addpath(genpath(fullfile(repo_root, 'src')));

% Construct the option struct.
% You can override parameters, e.g.:
%   opts = set_default_opts('N', 128, 'maxiter_opt', 40);
opts = set_default_opts('N', 64);

% Run the solver. The function will:
%   - run Algorithm 3 (outer loop + Manopt inner solver)
%   - save a .mat result file under results/ via result_filename(opts)
results = run_solver(opts);
