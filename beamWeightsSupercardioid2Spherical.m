function b_n = beamWeightsSupercardioid2Spherical(N)
%BEAMWEIGHTSSUPERCARDIOID2SPHERICAL Generate beamweights for supercardioids
%
%   For a specific order N, the supercardioid is the pattern that maximizes
%   the front-to-back power ratio. As no closed-form formula is known to the 
%   author that returns their coefficients directly in the SHD, beamweights
%   up to order 4 have been tabulated below. These have been converted from 
%   coefficients of differential form found in 
%
%       Elko, G.W., 2004. Differential microphone arrays. 
%       In Audio signal processing for next-generation multimedia communication systems (pp. 11-65). Springer.
%
%   Because the pattern is axisymmetric only the N+1 coefficients of m=0 
%   are returned.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSSUPERCARDIOID2SPHERICAL.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if N>4
    error('Coefficients available up to 4')
end

switch N
    case 1
        b_n = [1.2975    1.2975]';
    case 2
        b_n = [0.8372    0.9591    0.4680]';
    case 3
        b_n = [0.6267    0.7861    0.5021    0.1641]';
    case 4
        b_n = [0.5013    0.6674    0.4942    0.2314    0.0569]';
end
