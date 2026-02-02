function opts = set_default_opts(varargin)
%SET_DEFAULT_OPTS Centralized parameter setup for all scripts.
%
% Usage:
%   opts = set_default_opts('N', 128, 'maxiter_opt', 30);

opts = struct();

% ---- core problem size ----
opts.N   = 64;
opts.h   = 1.0;
opts.r_l = 100;

% ---- algorithm parameters ----
opts.mu          = 5;
opts.tau         = 2;
opts.sigma       = 1;
opts.maxiter     = 100;
opts.maxiter_opt = 30+10*(log2(opts.N)-6);

% ---- io ----
opts.root_dir    = fileparts(fileparts(mfilename('fullpath'))); % repo root
opts.data_dir    = fullfile(opts.root_dir, 'data');
opts.results_dir = fullfile(opts.root_dir, 'results');
opts.fix_file    = fullfile(opts.data_dir, 'fix.mat');

% ---- apply overrides (name-value pairs) ----
if mod(numel(varargin), 2) ~= 0
    error('set_default_opts:InvalidInput', ...
          'Overrides must be name-value pairs.');
end

for k = 1:2:numel(varargin)
    name = varargin{k};
    val  = varargin{k+1};

    if ~(ischar(name) || isstring(name))
        error('set_default_opts:InvalidName', ...
              'Override names must be char/string.');
    end
    name = char(name);
    opts.(name) = val;
end

% ---- dependent parameters ----
opts.m = log2(opts.N) - 3;

% ---- ensure folders exist (optional but helpful) ----
if ~exist(opts.data_dir, 'dir');    mkdir(opts.data_dir);    end
if ~exist(opts.results_dir, 'dir'); mkdir(opts.results_dir); end

end
