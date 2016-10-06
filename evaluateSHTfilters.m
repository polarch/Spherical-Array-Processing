function [cSH, lSH, WNG] = evaluateSHTfilters(M_mic2sh, H_array, fs, Y_grid, w_grid)
%EVALUATESHTFILTERS Evaluate frequency-dependent performance of SHT filters
%
%   The SHT filters can be evaluated in terms of how ideal are the SH
%   components that they generate. The evaluation here follows the metrics
%   introduced in 
%
%       Moreau, S., Daniel, J., Bertet, S., 2006, 
%       3D sound field recording with higher order ambisonics-objective measurements and validation of spherical microphone. 
%       In Audio Engineering Society Convention 120.
%
%   These are a) the spatial correlation coefficient between each ideal
%   spherical harmonic and the reconstructed pattern, evaluated at a dense
%   grid of directions, b) level difference between the mean spatial power
%   of the reconstructed pattern (diffuse power) over the one from an ideal
%   SH component. Ideally, correlaiton should be close to one, and the
%   level difference should be close to 0dB.
%
%   Additionally, the maximum amplification of all output SH components is
%   evaluated, through the maximum eigenvalue of the filtering matrix.
%
%   Inputs:
%       M_mic2sh:   (order+1)^2 x nMics x nBins SHT filtering matrix produced 
%           by one of the methods included in the library
%       H_array:    nBins x nMics x nGrid modeled or measured spherical array 
%           responses in a dense grid of nGrid directions
%       Y_grid:     nGrid x (order+1)^2 spherical harmonics matrix for the
%           nGrid directions of the evaluation grid
%       w_grid:     nGridx1 vector of integration weights for the grid 
%           points (optional, leave empty if not important)
%       fs:         sample rate for the output filters
%
%   Outputs:
%       H_filt: nSH x nMics x (Lfilt/2+1) returned filters in the frequency 
%           domain (half spectrum up to Nyquist)
%       h_filt: Lfilt x nMics x nSH impulse responses of the above filters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EVALUATESHTFILTERS.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sampling weights
nGrid = size(H_array, 3);
if (nargin<5) || isempty(w_grid)
    w_grid = 1/nGrid*ones(nGrid,1);
end

nBins = size(M_mic2sh,3);
nFFT = 2*(nBins-1);
f = (0:nFFT/2)*fs/nFFT;
order_sht = sqrt(size(Y_grid,2))-1;

% compute spatial correlations and integrated level difference between 
% ideal and reconstructed harmonics

% % power of ideal SHs over all grid directions (should be close to unity)
% for q = 1:(order_sht+1)^2
%     y_ideal_nm = Y_grid(:,q);
%     y_ideal_pwr(q) = (y_ideal_nm.*w_grid)'*y_ideal_nm;
% end
    
for kk=1:nBins
    H_kk = squeeze(H_array(kk,:,:));
    y_recon_kk = M_mic2sh(:,:,kk) * H_kk;
    for n=0:order_sht
        cSH_n = 0; % spatial correlation (mean per order)
        lSH_n = 0; % diffuse level difference (mean per order)
%         rSH_n = 0; % mean level difference (mean per order)
        for m = -n:n
            q = n^2+n+m+1;
            y_recon_nm = y_recon_kk(q,:).';
            y_ideal_nm = Y_grid(:,q);
            cSH_nm = (y_recon_nm.*w_grid)'*y_ideal_nm/ sqrt( (y_recon_nm.*w_grid)'*y_recon_nm );
            cSH_n = cSH_n + cSH_nm;
            lSH_nm = real( (y_recon_nm.*w_grid)'*y_recon_nm );
            lSH_n = lSH_n + lSH_nm;  
%             rSH_nm = sum((abs(y_recon_nm - y_ideal_nm).^2).*w_grid);
%             rSH_n = rSH_n + rSH_nm;
        end
        cSH(kk, n+1) = cSH_n/(2*n+1);
        lSH(kk, n+1) = lSH_n/(2*n+1);
%         rSH(kk, n+1) = rSH_n/(2*n+1);
    end
end

% maximum noise amplification of all filters in matrix
for kk=1:nBins
    eigM = eig(M_mic2sh(:,:,kk)'*M_mic2sh(:,:,kk));
    WNG(kk,1) = max(eigM);
end

% plots
for n=0:order_sht; str_legend{n+1} = num2str(n); end
figure
subplot(311)
semilogx(f, abs(cSH)), grid, legend(str_legend), axis([50 20000 0 1])
title('Spatial correlation')
subplot(312)
semilogx(f, 10*log10(lSH)), grid, legend(str_legend), axis([50 20000 -30 10])
title('Level difference')
subplot(313)
semilogx(f, 10*log10(WNG)), grid, set(gca, 'xlim',[50 20000])
title('Maximum amplification')
xlabel('Frequency (Hz)')
% subplot(414)
% semilogx(f, 10*log10(rSH)), grid, set(gca, 'xlim',[50 20000])
% title('MSE')
% xlabel('Frequency (Hz)')

end
