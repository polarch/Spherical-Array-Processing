function b_n = beamWeightsDifferential2Spherical(a_n)
%BEAMWEIGHTSDIFFERENTIAL2SPHERICAL Transform differential to spherical coefficients.
%
%   For a set of a_n coefficients describing a pattern of a differential 
%   array of the form
%   D(theta)=Sum_{n=1}^{N} a_n * cos^n(theta), 
%   generate the beamweights for the same pattern in the SHD. Since the 
%   pattern is axisymmetric only the N+1 coefficients of m=0 are returned.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSDIFFERENTIAL2SPHERICAL.M - 10/4/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% order N (or polynomial degree)
N = length(a_n)-1;

% build Legendre polynomial matrix
P_N = zeros(N+1);
for n = 0:N
    P_N(1:n+1,n+1) = returnLegePolyCoeffs(n);
end
% build the normalisation matrix
w = sqrt((2*(0:N)+1)/(4*pi));
W = diag(w);
W_inv = diag(1./w);
% convert coefficients
b_n = W_inv * inv(P_N) * a_n;
