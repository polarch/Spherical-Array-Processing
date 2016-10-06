function [P_music, est_dirs] = sphMUSIC(sphCOV, grid_dirs, nSrc)
%SPHMUSIC DoA estimation using MUSIC in the SHD
%   
%   This routine computes the MUSIC pseudo-spectrum directly in the SHD.
%   Subspace approaches such as MUSIC can offer higher spatial resolution 
%   than beamforming approaches such as the steered-response power, as long
%   as the source signals are not correlated between them and with the 
%   reverberant/diffuse sound.
%
%   Inputs:
%       sphCOV: (order+1)^2x(order+1)^2 covariance/correlation matrix of
%           SH signals
%       grid_dirs:  Kx2 directions of [azi elev] in rads that define a grid
%           for the power map. For easy plotting of the directional maps,
%           use grid2dirs.m to generate the directions
%       nSrc:   (optional) number of peaks to try to find, as potential DoA
%           estimates (Von-Mises peak-finding contributed by Dr. Sakari Tervo)
%
%   Outputs:
%       P_pwd:      Kx1 vector of output powers, evaluated at grid directions
%       est_dirs:   nSrcx2 [azi elev] of estimated directions from
%           peak-finding
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHMUSIC.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSH = size(sphCOV,1);
order = sqrt(nSH)-1;
% sorted eigenvalue decomposition
V = sorted_eig(sphCOV, 'descend');
% noise subspace
Vn = V(:,nSrc+1:end);

% get steering vectors for grid points
nGrid = size(grid_dirs,1);
grid_xyz = unitSph2cart(grid_dirs);
grid_dirs2 = [grid_dirs(:,1) pi/2-grid_dirs(:,2)]; % go from azi-elevation to azi-inclination
Y_grid = getSH(order, grid_dirs2, 'real');

% get MUSIC spectrum
P_music = zeros(nGrid,1);
for ng = 1:nGrid
    stVec = Y_grid(ng,:).';
    P_music(ng) = 1 / (stVec' * Vn * Vn' * stVec);
end

% peak finding, if asked
if nargout==2 && nargin==3
    
    kappa  = 50; % Von-Mises concentration factor
    P_minus_peak = P_music;
    est_dirs = zeros(nSrc, 2);
    for k = 1:nSrc
        [~, peak_idx] = max(P_minus_peak);
        est_dirs(k,:) = grid_dirs(peak_idx,:);
        VM_mean = grid_xyz(peak_idx,:); % orientation of VM distribution
        VM_mask = kappa/(2*pi*exp(kappa)-exp(-kappa)) * exp(kappa*grid_xyz*VM_mean'); % VM distribution
        VM_mask = 1./(0.00001+VM_mask); % inverse VM distribution
        P_minus_peak = P_minus_peak.*VM_mask;
    end
end

end
