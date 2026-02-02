function val=grad_S0(y, t, N, m, r_l, V, W, X0)
%GRAD_S0  Euclidean gradient for the initialization problem in init_S0.m
%
% This function provides the gradient of the following normalized least-squares
% objective used in init_S0:
%   f(y,t) = || X0 - H(y,t) ||^2 / ||X0||^2
%
% where H(y,t) is the hierarchical PSD representation (Eq. (25) in the paper),
%
% Intended usage (Manopt):
%   problemS.cost  = @(yt)  norm(X0 - H(yt.y,yt.t,...))^2 / norm(X0)^2;
%   problemS.egrad = @(yt)  grad_S0(yt.y,yt.t,...,V,W,X0);
%
% Inputs:
%   y   : complex matrix of size (3N) x (r_l*m)  
%   t   : complex vector of size (3N+1) x 1      
%   N   : system size 
%   m   : hierarchy layers
%   r_l : rank 
%   V,W : precomputed sparse operators 
%   X0  : target vectorized matrix (length (3N+1)^2), e.g. S0(:) in
%   init_S0.m
%
% Output:
%   val : struct with fields:
%         - val.y : Euclidean gradient w.r.t. y (size (3N) x (r_l*m))
%         - val.t : Euclidean gradient w.r.t. t (size (3N+1) x 1)

    dLdX=(X0-H(y,t,N,m,r_l))';
    
    r=reshape(dLdX*V,3*N,3*N*m);
    s=reshape(y.'*r,m*3*N*r_l*m,1);

    % Gradient w.r.t. y 
    dLdy=-2*reshape(W*s,3*N,r_l*m);

    % Gradient w.r.t. t:
    dLdt=-2*reshape(dLdX,3*N+1,3*N+1).'*t;
    
    val=struct('y',dLdy/norm(X0)^2,'t',dLdt/norm(X0)^2);
end
