function f = result_filename(opts)
%RESULT_FILENAME  Return the full path of the output result file for a given run.
%
% Usage:
%   f = result_filename(opts);
%
% Input:
%   opts : option struct (typically created by set_default_opts) that must contain:
%          - opts.results_dir : directory where result files are written (e.g., <repo>/results)
%          - opts.N           : problem size used to name the result file
%
% Output:
%   f : full file path to the result .mat file, formatted as:
%       <opts.results_dir>/<N>_benchmark.mat
%
% Example:
%   opts = set_default_opts('N', 64);
%   save(result_filename(opts), 'eta_P', 'eta_D', ...);

fname = sprintf('%d_benchmark.mat', opts.N);
f = fullfile(opts.results_dir, fname);
end

