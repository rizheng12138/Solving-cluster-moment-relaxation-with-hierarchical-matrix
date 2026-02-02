function val=H(y,t,N,m,r_l)
%H  Hierarchical representation of the PSD variable S (Eq. (25) in the paper).
%
% This function constructs a structured PSD matrix via a hierarchical 
% low-rank representation and returns it in vectorized form.
%
% Inputs:
%   y   : complex matrix of size (3N) x (r_l*m)
%         y is partitioned into m level-blocks y_i = y(:, r_l*(i-1)+1 : r_l*i),
%         each of size (3N) x r_l.
%   t   : complex vector of size (3N+1) x 1 
%   N   : system size 
%   m   : hierarchy layers 
%   r_l : rank parameter per level
%
% Output:
%   val : vectorized matrix H (length (3N+1)^2).
%
% Notes:
%   - Each term v*v' is PSD, and t*t' is PSD, hence the resulting S is PSD.

    val=zeros(3*N+1,3*N+1);

    % H^{(2)}(y)
    for i=1:1:m
        len=3*N/2^(i-1);

        y_i=y(:,r_l*(i-1)+1:r_l*i);

        for j=1:1:2^(i-1)
            v=y_i((j-1)*len+1:j*len,:);

            val((j-1)*len+1:j*len,(j-1)*len+1:j*len)= ...
                val((j-1)*len+1:j*len,(j-1)*len+1:j*len)+v*v';
        end
    end

    % Add the final rank-1 PSD term t*t' 
    val=val+t*t';

    val=val(:);
end
