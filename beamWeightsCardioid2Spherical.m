function b_n = beamWeightsCardioid2Spherical(N)
%BEAMWEIGHTSCARDIOID2SPHERICAL Generate spherical coefficients for cardioids
%
%   For a specific order N of a higher order cardioid of the form
%   D(theta)=(1/2)^N * (1+cos(theta))^N, generate the beamweights for the
%   same pattern in the SHD. Because the pattern is axisymmetric only the N+1
%   coefficients of m=0 are returned.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSCARDIOID2SPHERICAL.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The coefficients can be derived by the binomial expansion of the
% cardioid function
b_n = zeros(N+1,1);
for n=0:N
    b_n(n+1) = sqrt(4*pi*(2*n+1)) * factorial(N)*factorial(N+1)/(factorial(N+n+1)*factorial(N-n))/(N+1);
end
