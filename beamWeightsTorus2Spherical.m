function b_n = beamWeightsTorus2Spherical(N)
%BEAMWEIGHTSTORUS2SPHERICAL Generate beamweights for a raised torus
%
%   The tabulated beamweights below describe toroidal patterns of the form
%   D(theta) = |sin(\theta)|^N, using up to 4th-order components.
%   Note that N=2 and N=4 are exactly band-limited to 2nd and 4th-order
%   respectively, while the N=1 and N=3 extend to higher-orders with
%   decaying energy, and hence the weights below are just an approximation.
%   Because the pattern is axisymmetric only the N+1 coefficients of m=0 
%   are returned.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSTORUS2SPHERICAL.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch N
    case 1
        b_n = [2.7842    0.0000   -0.7781    0.0000   -0.1303];
    case 2
        b_n = [2.3633   -0.0000   -1.0569];
    case 3
        b_n = [2.0881   -0.0000   -1.1673    0.0000    0.1468];
    case 4
        b_n = [1.8906   -0.0000   -1.2079    0.0000    0.2701];
end
