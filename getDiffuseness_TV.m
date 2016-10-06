function diff = getDiffuseness_TV(i_vecs)
%GETDIFFUSENESS_TV Compute diffuseness from temporal variation of intensity
%vectors
%   
%   This routine estimates diffuseness using the temporal variation of 
%   intensity vectors, as proposed in 
%       
%       Ahonen, J. and Pulkki, V., 2009. 
%       Diffuseness estimation using temporal variation of intensity vectors. 
%       In 2009 IEEE Workshop on Applications of Signal Processing to Audio and Acoustics (WASPAA).
%
%   also termed coefficient of variation in other works. It is based on the
%   relation sqrt( 1 - ||<Ia>||/<||Ia||> ), where Ia is the active intensity 
%   vector. It is similar to the spherical variance estimator, but it takes 
%   into account the magnitude of the DoA vectors, hence giving more weight 
%   to strong sound estimates instead of low noisy ones. It is also flexible 
%   in the sense that alternative DoA vectors with magnitude that relates 
%   to the energy of the sound-field can be used instead of the intensity
%   vector.
%
%   Inputs:
%       i_vecs:  Kx3 matrix of DoA vectors, for K temporal observations
%
%   Outputs:
%       diff:       diffuseness value
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GETDIFFUSENESS_TV.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mean of active intensity vectors
mean_i_vecs = mean(i_vecs, 1);
% norm of mean vector
norm_mean_i_vecs = sqrt(sum(mean_i_vecs.^2, 2));

% norm of intensity vectors
norm_i_vecs = sqrt(sum(i_vecs.^2, 2));
% mean of norms
mean_norm_i_vecs = mean(norm_i_vecs);

% diffuseness
diff = sqrt( 1 - norm_mean_i_vecs/mean_norm_i_vecs );

end
