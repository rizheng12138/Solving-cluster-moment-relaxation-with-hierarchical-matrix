function f = data_filename(opts)
%DATA_FILENAME  Return the full path of the prepared data file for a given run.
%
% Usage:
%   f = data_filename(opts);
%
% Input:
%   opts : option struct (typically created by set_default_opts) that must contain:
%          - opts.data_dir : directory where data files are stored (e.g., <repo>/data)
%          - opts.N        : problem size used to name the dataset
%
% Output:
%   f : full file path to the dataset corresponding to opts.N, formatted as:
%       <opts.data_dir>/<N>_data.mat
%
% Example:
%   opts = set_default_opts('N', 64);
%   load(data_filename(opts));   % loads data/64_data.mat

f = fullfile(opts.data_dir, sprintf('%d_data.mat', opts.N));
end

