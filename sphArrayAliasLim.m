function [f_alias, cond_N] = sphArrayAliasLim(R, Nmic, maxN, mic_dirs_rad, mic_weights)
%SPHARRAYALIASLIM Get estimates of the aliasing limit of a spherical array
%   
%   First estimate takes into account only the radius and a nominal order
%   that the array is expected to support, it is the simplest one and it
%   expresses the kR = maxN rule.
%   The second estimate is based on the number of microphones, and it can
%   be more relaxed than the first, if the nominal supported order is less
%   than maxN<floor(sqrt(Nmic)-1).
%   The third takes into account microphone numbers and directions, and it
%   is based on the orthogonality of the SH matrix for the microphone
%   positions expressed through the condition number.
%
%   Inputs:
%       R:      radius of spherical array
%       Nmic:   number of microphones
%       maxN:   expected maximum supported order of the array
%       mic_dirs_rad:   nMicx2 matrix of [azi elev] directions of the
%           microphones on the array
%       mic_weights: (optional) nMicx1 vector of weights used to improve
%           orthogonality of the SH transform
%
%   Outputs:
%       f_alias:    the alising limit estimates
%       cond_N:     the condition number of all orders up to one above 
%           floor(sqrt(Nmic)-1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHARRAYALIASLIM.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 343;

% conventional kR=N assumption
f_alias(1) = c*maxN/(2*pi*R);

% based on the floor of the number of microphones, uniform arrangement
f_alias(2) = c*floor(sqrt(Nmic)-1)/(2*pi*R);

% based on condition number of the SHT matrix
aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

maxOrder = ceil(sqrt(Nmic)-1);
if (nargin<5), mic_weights = []; end

cond_N = checkCondNumberSHT(maxOrder, aziElev2aziIncl(mic_dirs_rad), 'real', mic_weights);
% plot(cond_N)
trueMaxOrder = find(cond_N<10^4, 1, 'last')-1;
f_alias(3) = c*trueMaxOrder/(2*pi*R);

end
