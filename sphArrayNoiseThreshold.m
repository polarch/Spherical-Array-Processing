function f_lim = sphArrayNoiseThreshold(R, Nmic, maxG_db, maxN, arrayType, dirCoeff)
%SPHARRAYNOISETHRESHOLD Returns freq.limits for noise amplification in an SMA
%
%   SPHARRAYNOISETHRESHOLD returns the frequencies that the noise in the 
%   output channels of a SMA, after performing the SHT and equalization of 
%   the output signals, reaches a certain user-defined threshold maxG_db.
%   The frequencies are computed only at the lower range of each order, 
%   where its response decays rapidly, ignoring for example the nulls of an 
%   open array at the higher frequencies. The estimation of the limits are
%   based on a linear approximation of the log-log response found e.g. in 
%
%       Sector-based Parametric Sound Field Reproduction in the Spherical Harmonic Domain
%       A Politis, J Vilkamo, V Pulkki
%       IEEE Journal of Selected Topics in Signal Processing 9 (5), 852 - 866
%
%   Inputs:
%         R:          array radius
%         Nmic:       number of microphones
%         maxG_db:    max allowed amplification for the noise level, 
%               maxG_db = 20*log10(maxG)
%         maxN:       maximum order that the array supports
%         arrayType:  'open', 'rigid' or 'directional' for an open or rigid spherical 
%               array of omnidirectional microphones, or for an array of
%               first-order directional microphones
%         dirCoeff: (optional) if arrayType=='directional', then dirCoeff 
%               expresses the first-order directivity of microphones
%               from 0-1 (1:omni, 0.5:cardioid, 0:dipole)
% 
%         f_lim:    frequency points that the threhsold is reached, 
%               for orders 1:maxN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHARRAYNOISETHRESHOLD.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    dirCoeff = [];
end

c = 343;
f_lim = zeros(1,maxN);
for n=1:maxN
    bn = sphModalCoeffs(n, 1, arrayType, dirCoeff)/(4*pi);
    bn = bn(end);
    maxG = 10^(maxG_db/10);
    kR_lim = (maxG*Nmic*abs(bn).^2)^(-10*log10(2)/(6*n));
    f_lim(n) = kR_lim*c/(2*pi*R);
end

end
