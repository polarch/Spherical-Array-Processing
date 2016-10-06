function W_nullpwd = sphNullformer_pwd(order, beam_dirs)
%SPHNULLFORMER_PWD Beamweights for a PWD beamformer at a specified direction, 
%with nulls at others
%   
%   For a set of K directions, computes the beamforming weights that creates
%   a PWD beamformer at each of the direction in the set, while placing
%   nulls at the rest.
%
%   Inputs:
%       order:  order of SH signals
%       beam_dirs:  Kx2 [azi elev] directions
%
%   Outputs:
%       W_nullpwd:  (order+1)^2xK matrix of beamweights, one set per
%           column for each beamforming direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHNULLFORMER_PWD.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    beam_dirs2 = [beam_dirs(:,1) pi/2-beam_dirs(:,2)]; % convert from azi-elev to azi-incl

    % steering vectors in the SHD
    Y_nullpwd = getSH(order, beam_dirs2, 'real');

    % beamform at each direction and nullform at the rest
    W_nullpwd = pinv(Y_nullpwd);
end
