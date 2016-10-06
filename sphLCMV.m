function w_lcmv = sphLCMV(sphCOV, constraint_dirs, constraints)
%SPHLCMV Get beamforming weights for a LCMV beamformer in the SHD
%   
%   Computes the beamweights for a linearly constrained minimum variance
%   (LCMV) beamformer using SH signals.
%
%   Inputs:
%       sphCOV:     (order+1)^2x(order+1)^2 covariance/correlation matrix
%           of SH signals
%       constraint_dirs:    Kx2 [azi elev] directions that the beamformer
%           should achieve the user-specified constraints
%       constraints:        Kx1 constraint values
%
%   Outputs:
%       w_lcmv: (order+1)^2x1 vector of beamforming weights
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHLCMV.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSH = size(sphCOV,1);
order = sqrt(nSH)-1;

% steering vectors for distortionless constraints
constraint_dirs2 = [constraint_dirs(:,1) pi/2-constraint_dirs(:,2)]; % from azi-elevation to azi-inclination
Y_lcmv = getSH(order, constraint_dirs2, 'real');

% LCMV weights vector
stVecMtx = Y_lcmv.';
invA_B = sphCOV\stVecMtx;
B_invA_B = stVecMtx' * invA_B;
w_lcmv = (invA_B / B_invA_B) * constraints;

end
