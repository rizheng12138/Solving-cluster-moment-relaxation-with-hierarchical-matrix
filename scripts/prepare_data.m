% prepare_data.m
% -------------------------------------------------------------------------
% Prepare and save problem data for the 1D TFI benchmark instance.

clear; clc;
this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(this_file));  
addpath(genpath(fullfile(repo_root, 'src')));

% Set experiment options. Default here uses N=64 for quick initialization.
% You can override settings by adding name-value pairs, e.g.:
%   opts = set_default_opts('N', 128, 'r_l', 120);
opts = set_default_opts('N', 64);

% Unpack frequently used parameters for readability.
N=opts.N;
h=opts.h;
r_l=opts.r_l;
m=opts.m;
mu=opts.mu;
tau=opts.tau;
maxiter=opts.maxiter;
maxiter_opt=opts.maxiter_opt;
sigma=opts.sigma;

% Compute true ground-state energy
if h==1.0
    true_E0=1.2732*N;
elseif h==0.5
    true_E0=1.0635*N;
elseif h==1.5
    true_E0=1.6719*N;
end

% Compute J in primal objective tr(J*M)
J=get_J(N,h);

% Compute projection operator R needed for run_solver.m
T=sparse((3*N+1)^2,(3*N+1)^2,1,(3*N+1)^2,(3*N+1)^2);
u=sparse(ones(1,3*N+1),(0:3*N)*(3*N+1)+(1:3*N+1),-1,1,(3*N+1)^2);

row1=zeros(1,15*N);
col1=zeros(1,15*N);

temp_row2=zeros(1,12*N);
temp_col2=zeros(1,12*N);

for i=1:N
    row1(15*i-14:15*i)=[
        (3*i-3)*(3*N+1)+3*i-2,(3*i-3)*(3*N+1)+3*i-1,(3*i-3)*(3*N+1)+3*i, ...
        (3*i-2)*(3*N+2),(3*i-1)*(3*N+1)+(3*i-2),(3*i-3)*(3*N+1)+(3*i-1),...
        (3*i-2)*(3*N+1)+3*i-2,(3*i-2)*(3*N+1)+3*i-1,(3*i-2)*(3*N+1)+3*i, ...
        (3*i-1)*(3*N+2),(3*i-3)*(3*N+1)+3*i,(3*i-2)*(3*N+1)+3*i,...
        (3*i-1)*(3*N+1)+3*i-2,(3*i-1)*(3*N+1)+3*i-1,(3*i-1)*(3*N+1)+3*i];
    col1(15*i-14:15*i)=[
        (3*i-3)*(3*N+1)+3*i-2,(3*i-3)*(3*N+1)+3*i-1,(3*i-3)*(3*N+1)+3*i, ...
        (3*i-3)*(3*N+1)+(3*i-1),(3*i-3)*(3*N+1)+3*i,(3*i-2)*(3*N+2),...
        (3*i-2)*(3*N+1)+3*i-2,(3*i-2)*(3*N+1)+3*i-1,(3*i-2)*(3*N+1)+3*i, ...
        (3*i-2)*(3*N+1)+3*i,(3*i-1)*(3*N+1)+3*i-2,(3*i-1)*(3*N+2),...
        (3*i-1)*(3*N+1)+3*i-2,(3*i-1)*(3*N+1)+3*i-1,(3*i-1)*(3*N+1)+3*i];

    temp_row2(12*i-11:12*i)=repmat(([1 2 0 2 0 1]+3*(i-1))*(3*N+1)+([3 2 3 1 2 1]+3*(i-1)),1,2);
    temp_col2(12*i-11:12*i)=[repelem((3*i-2:3*i)*(3*N+1),2),repelem(3*N*(3*N+1)+3*i-2:3*N*(3*N+1)+3*i,2)];
end

row2=[temp_row2,temp_col2];
col2=[temp_col2,temp_row2];
val1=repmat([1,0.75,0.75,0.25,0.25,0.25,0.75,1,0.75,0.25,0.25,0.25,0.75,0.75,1],1,N);
temp_val2=repmat([1,-1,-1,1,1,-1,1,-1,-1,1,1,-1],1,N);
val2=[0.25i*temp_val2,-0.25i*temp_val2];

row3=[(3*N+1)*3*N+1:(3*N+2)*3*N,3*N+1:3*N+1:3*N*(3*N+1),(3*N+1)*3*N+1:(3*N+2)*3*N,(3*N+1)*(1:3*N)];
col3=[(3*N+1)*3*N+1:(3*N+2)*3*N,3*N+1:3*N+1:3*N*(3*N+1),(3*N+1)*(1:3*N),(3*N+1)*3*N+1:(3*N+2)*3*N];
val3=[0.75*ones(1,6*N),-0.25*ones(1,6*N)];

T=T+sparse([row1,row2,row3],[col1,col2,col3],[val1,val2,val3],(3*N+1)^2,(3*N+1)^2);

clear col1 col2 col3 row1 row2 row3 temp_row2 temp_col2 temp_val2 val1 val2 val3;

[A,B] = ndgrid(1:3*N+1, 1:3*N+1);
rows = (A(:)-1) * (3*N+1) + B(:);
cols = (B(:)-1) * (3*N+1) + A(:);
K = sparse(rows, cols, 1, (3*N+1)^2, (3*N+1)^2);

index = reshape(((3*(1:N) - [3; 2; 1]) * (3*N + 1) + 3*(1:N) - 3) + reshape(1:3, 1, 1, 3),1,9*N);
index=setdiff(1:(3*N+1)*3*N,cat(2,3*N+1:3*N+1:(3*N+1)*3*N,index));
T=T+0.5*(speye((3*N+1)^2,(3*N+1)^2)-K)*sparse(index,index,1,(3*N+1)^2,(3*N+1)^2);

clear A B rows cols K index;

R=speye((3*N+1)^2,(3*N+1)^2)-T;
clear T;
clear K L l i k row_idx col_idx idx offsets OffsetGrid IdxGrid lin_idx;

% Precompute constants for feasibility measures and dual objective
% cD: scaling constant for dual feasibility measure
cD=1/(1+norm(J));

% idx_u: linear indices selecting diagonal entries of a (3N+1)x(3N+1) matrix
% (used to form dual_obj = sum(C(idx_u)))
idx_u=(0:3*N)*(3*N+1)+(1:3*N+1);

% eigopts: options for eigs; v0 will be used later for warm-start
eigopts.isreal = true;


% Allocate history arrays (outer iterations)
eta_P=zeros(maxiter,1);   % primal feasibility measure history
eta_D=zeros(maxiter,1);   % dual feasibility measure history
eta_g=zeros(maxiter,1);   % duality gap measure history
energy=zeros(maxiter,1);  % primal objective history


% Initialize primal variable M
% M represents a (3N+1)x(3N+1) matrix but is stored as a vector for efficiency.
M = eye(3*N+1);
M=M(:);

% save the pre-generated data
save(data_filename(opts), '-v7.3');
fprintf('Complete!\n');

