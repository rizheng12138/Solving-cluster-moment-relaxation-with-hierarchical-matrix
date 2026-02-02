function f = s0_filename(opts)
%S0_FILENAME  Return the full path of the saved initialization (S0) file.
%
% Usage:
%   f = s0_filename(opts);
%
% Input:
%   opts : option struct (typically created by set_default_opts) that must contain:
%          - opts.data_dir : directory where data/initialization files are stored
%          - opts.N        : problem size used to name the initialization file
%
% Output:
%   f : full file path to the initialization .mat file, formatted as:
%       <opts.data_dir>/<N>_S0.mat
%
% Example:
%   opts = set_default_opts('N', 64);
%   load(s0_filename(opts));   % loads data/64_S0.mat

f = fullfile(opts.data_dir, sprintf('%d_S0.mat', opts.N));
end

