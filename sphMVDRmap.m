function [P_mvdr, est_dirs] = sphMVDRmap(sphCOV, grid_dirs, nSrc)
%SPHMVDRMAP Steered-response power map using a MVDR beamformer
%   
%   This routine steeres a MVDR beamformer at K grid directions
%   and computes the power output. This power map can be used for DoA
%   estimation. The rejection of diffuse sound and interferes of MVDR
%   results in higher spatial resolution on the power map than with
%   signal-independent beamforming such as a regular (PWD) one.
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
% SPHMVDRMAP.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MVDR power map
nGrid = size(grid_dirs,1);
grid_xyz = unitSph2cart(grid_dirs);
P_mvdr = zeros(nGrid,1);
for ng = 1:nGrid
    w_mvdr = sphMVDR(sphCOV, grid_dirs(ng,:));
    P_mvdr(ng) = w_mvdr' * sphCOV * w_mvdr;
end

% peak finding, if asked
if nargout==2 && nargin==3
    
    kappa  = 20; % Von-Mises concentration factor
    P_minus_peak = P_mvdr;
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
