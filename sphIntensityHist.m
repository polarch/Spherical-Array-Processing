function [I_hist, est_dirs] = sphIntensityHist(i_xyz, grid_dirs, nSrc)
%SPHINTENSITYHIST Form acoustic intensity histograms and estimate DoAs
%   
%   This routine takes a number of acoustic intensity measurements, sampled
%   across time, or frequency bands, and computes histograms of their DoAs,
%   weighted by the magnitudes of the vectors.
%
%   Inputs:
%       i_xyz:  Kx3 matrix of [i_x i_y i_z] values for K intensity
%           observations
%       grid_dirs:  Kx2 directions of [azi elev] in rads that define a grid
%           for the histograms. For easy plotting of the direcitonal maps,
%           use grid2dirs.m to generate the directions
%       nSrc:   (optional) number of peaks to try to find, as potential DoA
%           estimates (Von-Mises peak-finding contributed by Dr. Sakari Tervo)
%
%   Outputs:
%       I_hist:     Kx1 histogram, evaluated at grid directions
%       est_dirs:   nSrcx2 [azi elev] of estimated directions from
%           peak-finding
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHINTENSITYHIST.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nGrid = size(grid_dirs,1);
nSamples = size(i_xyz,1);
% find closest direction on grid for each sample
grid_xyz = unitSph2cart(grid_dirs);
i_dirs_idx = zeros(nSamples,1);
for ns = 1:nSamples
    [~, i_dirs_idx(ns)] = min( sum( (ones(nGrid,1)*i_xyz(ns,:) -  grid_xyz).^2, 2) );
end
% take intensity magnitude as histogram weighting
I_mag = sum( i_xyz.^2, 2);
I_hist = zeros(nGrid,1);
for ng=1:nGrid
    I_hist(ng) = sum(I_mag(i_dirs_idx==ng));
end

% peak finding, if asked
if nargout==2 
    if nargin<3, nSrc = 1; end
    
    kappa  = 20; % Von-Mises concentration factor
    P_minus_peak = I_hist;
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
