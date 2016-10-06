function W_mvdr = sphMVDR(sphCOV, beam_dirs)
%SPHMVDR Get beamforming weights for an MVDR beamformer in the SHD
%   
%   Computes the beamweights for a minimum variance distortionless response
%   (MVDR) beamformer using SH signals. If multiple directions are passed
%   for distortionless response one set of weights is computed for each one
%   of them.
%
%   Inputs:
%       sphCOV:     (order+1)^2x(order+1)^2 covariance/correlation matrix
%           of SH signals
%       beam_dirs:  Kx2 [azi elev] beamforming directions for
%           distortioneless constraint of MVDR
%
%   Outputs:
%       W_mvdr:     (order+1)^2xK matrix of beamweights, one set per
%           column for each beamforming direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHMVDR.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSH = size(sphCOV,1);
order = sqrt(nSH)-1;
nBeams = size(beam_dirs,1);

% MVDR weights matrix
W_mvdr = zeros(nSH, nBeams);
% steering vectors for distortionless constraints
beam_dirs2 = [beam_dirs(:,1) pi/2-beam_dirs(:,2)]; % from azi-elevation to azi-inclination
Y_mvdr = getSH(order, beam_dirs2, 'real');

for nb=1:nBeams
    stVec = Y_mvdr(nb,:).';
    invA_b = sphCOV\stVec;
    b_invA_b = stVec' * invA_b;
    W_mvdr(:,nb) = invA_b / b_invA_b;
end

end
