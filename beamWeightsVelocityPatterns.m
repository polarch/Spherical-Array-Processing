function velCoeffs = beamWeightsVelocityPatterns(b_n, orientation, A_xyz, basisType)
%BEAMWEIGHTSVELOCITYPATTERNS Generate weighted velocity patterns for a
%given sector pattern
%
%   If the sound-field is weighted with an axisymmetric spatial distribution 
%   described by the N+1 SH coefficients b_n, then the beamweights
%   capturing the velocity signals for the weighted sound-field are of
%   an order one higher than the weighting pattern, and can be derived from
%   it. This type of beamforming has some applications on spatial sound
%   reproduction and acoustic analysis, see
%
%       Politis, A. and Pulkki, V., 2016. 
%       Acoustic intensity, energy-density and diffuseness estimation in a directionally-constrained region. 
%       arXiv preprint arXiv:1609.03409.
%
%   and
%
%       Politis, A., Vilkamo, J. and Pulkki, V., 2015. 
%       Sector-Based Parametric Sound Field Reproduction in the Spherical Harmonic Domain. 
%       IEEE Journal of Selected Topics in Signal Processing, 9(5), pp.852-866.
%
%   Inputs:
%       b_n:    N+1 beamweights in the SHD for axisymmetric weighting
%           sector pattern
%       orientation:    [azi elev] of look-up direction of sector pattern
%       A_xyz:  pre-computed matrices for the conversion (see
%           computeVelCoeffsMtx.m)
%       basisType:  'real or 'complex' to get the beamweights for the
%           respective SH signal convention
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSVELOCITYPATTERNS.M - 7/06/2015
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(b_n)-1;

azi = orientation(1);
polar = pi/2-orientation(2);
c_nm = rotateAxisCoeffs(b_n, polar, azi, 'complex');

% get the velocity coeffs
x_nm = A_xyz(1:(N+2)^2, 1:(N+1)^2, 1)*c_nm;
y_nm = A_xyz(1:(N+2)^2, 1:(N+1)^2, 2)*c_nm;
z_nm = A_xyz(1:(N+2)^2, 1:(N+1)^2, 3)*c_nm;

switch basisType
    case 'complex'
        velCoeffs = [x_nm y_nm z_nm];
    case 'real'
        % convert to real SH coefficients to use with the real signals
        velCoeffs = real(complex2realCoeffs([x_nm y_nm z_nm]));
end
