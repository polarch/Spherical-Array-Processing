function [V_sorted, S_sorted] = sorted_eig(X, direction)
%SORTED_EIG Compute the sorted eigenvalue decomposition of a square matrix
%
%   Inputs:
%       X:          matrix to decompose
%       direction:  'ascend','descend' for ordering of eigenvectors and
%           eigenvalues
%
%   Outputs:
%       V_sorted:   matrix of eigenvectors sorted according to the
%           eigenvalues
%       S_sorted:   diagonal matrix of sorted eigenvalues
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SORTED_EIG.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   

if nargin<2, direction = 'descend'; end

if size(X,1) ~= size(X,2)
    error('input matrix should be square')
end

[V, S] = eig(X);
[~, perm] = sort(diag(S), 1, direction);
S_sorted = S(perm, perm); V_sorted = V(:, perm);

end
