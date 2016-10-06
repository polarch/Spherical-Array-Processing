function [g2, g2_lin] = sphArrayNoise(R, Nmic, maxN, arrayType, f)
%SPHARRAYNOISE Returns noise amplification curves of a spherical mic. array
%
%   SPHARRAYNOISE returns the noise amplification curves with frequency and
%   order of the output SH components, for a theoretical microphone array
%   with assumed uniformly distributed microphones. For an assumed unity
%   noise variance of the microphones, these curves would express the power
%   spectral densities of the noise at the output channels after the SHT
%   and perfect equalization of the components.
%
%   R:          array radius
%   Nmic:       number of microphones
%   maxN:       maximum order that the array supports
%   arrayType:  'open' or 'rigid' for an open sphere array or for
%               microphones mounted on a rigid baffle
%   f:          frequency points to evaluate the curve
%   
%   g2:         the noise amplification curve
%   g2_lin:     an approximation of the curves at low frequencies showing 
%               the linear behaviour in the log-log axis, with a 6n dB
%               slope
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHARRAYNOISE.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% frequency axis
c = 343;
kR = 2*pi*f*R/c;
% modal responses
bN = sphModalCoeffs(maxN, kR, arrayType)'/(4*pi);
% noise power response
g2 = 1./(Nmic*abs(bN).^2).';

% approximate linearly
p = -(6/10)*(1:maxN)/log10(2);
bN_lim0 = sphModalCoeffs(maxN, 1, arrayType)/(4*pi);
a = 1./(Nmic*abs(bN_lim0(2:end)).^2);

g2_lin = zeros(length(kR), maxN);
for n=1:maxN
    g2_lin(:,n) = a(n)*kR.^p(n);
end

end
