function d_n = beamWeightsMaxEV(order)
%BEAMWEIGHTSMAXEV Generate spherical beamweights for maximum energy-vector patterns
%
%   Generate the beamweights for the a maximum energy-vector beampattern in
%   the SHD. This pattern originates from ambisonic-related research and it
%   maximizes the ambisonic energy-vector, which is essentially the
%   directional centroid of the squared pattern. IT can also be seen as the
%   pattern that maximizes the acoustic intensity vector of a diffuse
%   field weighted with this pattern. In practice it is almost the same as
%   a supercardioid that maximizes front-back power ratio for a certain order, 
%   and it can be used as such. Because the pattern is axisymmetric only the 
%   N+1 coefficients of m=0 are returned. Details for their theory can be
%   found e.g. in 
%
%       Zotter, F., Pomberger, H. and Noisternig, M., 2012. 
%       Energy-preserving ambisonic decoding. Acta Acustica united with Acustica, 98(1), pp.37-47.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSMAXEV.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_n = zeros(order+1,1);
for n=0:order
    temp = legendre(n, cos(2.4068/(order+1.51)));
    d_n(n+1) = sqrt((2*n+1)/(4*pi))*temp(1);
end
% normalize to unity response on look-direction
norm = sum(d_n.*sqrt((2*(0:order)'+1)/(4*pi)));
d_n = d_n/norm;
