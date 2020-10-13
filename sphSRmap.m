function [P_sr, est_dirs] = sphSRmap(shsig, p, A_grid, regValue, stopValue, maxIter, grid_dirs, nSrc)
%SPHSRMAP Spatial power map based on sparse recovery
%   
%   This routine performs sparse recovery using the iterative reweighted
%   least squares (IRLS) algorithm, with an L_{p,2} norm, on either time-domain
%   SH signals, or frquency-domain SH signals (at a single frequency) using 
%   a few consecutive STFT frames.
%
%   Inputs:
%       shsig:  SH signals of order N, with channel dimension first
%       p:  sparsity exponent p<=1, of the L_{p,2} norm
%       A:  overcomplete dictionary (in CS jargon) of steering vectors for
%       a dense grid of K directions, with A of size MxK.
%       nSrc:   (optional) number of peaks to try to find, as potential DoA
%           estimates (Von-Mises peak-finding contributed by Dr. Sakari Tervo)
%
%   Outputs:
%       P_sr:       Kx1 vector of output powers, evaluated at grid directions
%       est_dirs:   nSrcx2 [azi elev] of estimated directions from
%           peak-finding
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHSRMAP.M - 13/5/2019
% Archontis Politis, archontis.politis@tuni.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid_xyz = unitSph2cart(grid_dirs);

% compute solution (get only the resulting powers of the solution for each
% grid point)
[~,~,P_sr] = sparse_solver_irls(p, A_grid, shsig, regValue, stopValue, maxIter);

% peak finding, if asked
if nargout==2
    
    kappa  = 20; % Von-Mises concentration factor
    P_minus_peak = P_sr;
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