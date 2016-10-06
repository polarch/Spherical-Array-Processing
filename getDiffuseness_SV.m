function diff = getDiffuseness_SV(i_vecs)
%GETDIFFUSENESS_SV Compute the DoA spherical variance diffuseness 
%   
%   This routine estimates diffuseness based on the assumption that DoA 
%   estimates should have a small variance around the source DoA, if the
%   direct-to-diffuse ratio is high, and be completely uniformly
%   distributed if the sound is fully diffuse. This variation can be
%   captured on the spherical variance of the DoA estimates. The method can
%   utilize intensity vectors, or any other vectors that indicate DoAs. The
%   statistics can be computed from multiple temporal observations or
%   across different frequency bands. This estimator has been used e.g. in
%
%       Politis, A., Delikaris-Manias, S. and Pulkki, V., 2015. 
%       Direction-of-arrival and diffuseness estimation above spatial aliasing for symmetrical directional microphone arrays. 
%       In 2015 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP).
%
%   Inputs:
%       i_vecs:  Kx3 matrix of DoA vectors, for K temporal or frequency
%           observations
%
%   Outputs:
%       diff:       diffuseness value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GETDIFFUSENESS_SV.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unit DoA vectors
doa_vecs = i_vecs ./ (sqrt(sum(i_vecs.^2,2))*ones(1,3));
% mean DoA vector
mean_doa_vec = mean(doa_vecs,1);
% spherical variance
diff = 1 - sqrt(sum(mean_doa_vec.^2, 2));

end
