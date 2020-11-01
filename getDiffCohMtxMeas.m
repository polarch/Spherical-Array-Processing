function M_diffcoh = getDiffCohMtxMeas(H_array, w_grid)
%GETDIFFCOHMTXMEAS Returns the diffuse coherence matrix of the array
%   GETDIFFCOHMTXMEAS takes the complex measured array response for a grid
%   of directions around the array, and computes its Spatial Covariance
%   Matrix for a theoretical isotropic diffuse field of unit power. We
%   refer to this as the diffuse field coherence matrix D of the array where
%   each element (i,j) of the matrix represents the complex coherence
%   between channel i and j, for diffuse excitation. Technically, the
%   element should be also normalized with the signal powers to correspond 
%   to the standard coherence formula, as: D_ij/sqrt(DII*Djj). This matrix
%   is very useful for diffuse field modeling and synthesis, diffuse0-field
%   equalization techniques, and spatial filtering that involves models of
%   ambient noise as diffuse, for suppresion.
%
%   Inputs:
%       H_array:    nBins x nMics x nGrid modeled or measured spherical array 
%           responses in a dense grid of nGrid directions
%       w_grid:     nGridx1 vector of integration weights for the grid 
%           points (optional, leave empty if the points are close to uniform)
%
%   Outputs:
%       M_diffcoh: nMics x nMics x nBins diffuse-field coherence matrix of
%           the array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GETDIFFCOHMTXMEAS.M - 8/6/2017
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nMics = size(H_array,2);
nBins = size(H_array,1);
nGrid = size(H_array,3);

M_diffcoh = zeros(nMics,nMics,nBins);
if (nargin<2) || isempty(w_grid)
    w_grid = ones(nGrid,1)/nGrid;
end
W = diag(w_grid);
for nb=1:nBins
   H_nb = permute(H_array(nb,:,:),[2 3 1]);
   M_diffcoh(:,:,nb) = H_nb*W*H_nb';
end

end
