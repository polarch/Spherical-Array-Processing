function b_n = beamWeightsHypercardioid2Spherical(N)
%BEAMWEIGHTSHYPERCARDIOID2SPHERICAL Generate spherical beamweights for hypercardioids
%
%   The hypercardioid is the pattern that maximizes the directivity-factor
%   for a certain SH order N. The hypercardioid is also the plane-wave
%   decomposition beamformer in the SHD, also called 'regular' because the
%   beamweights are just the SH values on the beam-direction. 
%   Since the pattern is axisymmetric only the N+1 coefficients of m=0 are 
%   returned.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSHYPERCARDIOID2SPHERICAL.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_n = 4*pi/(N+1)^2 * getSH(N,[0 0],'real');
b_n = zeros(N+1,1);
for n=0:N
    b_n(n+1) = c_n((n+1)^2 -n);
end
