function [H_filt, h_filt] = arraySHTfiltersMeas_regLSHD(H_array, order_sht, grid_dirs_rad, w_grid, nFFT, amp_threshold)
%ARRAYSHTFILTERSMEAS_REGLSHD Generate SHT filters based on measured responses 
%(regularized least-squares in the SHD)
%
%   Generate the filters to convert microphone signals from a spherical
%   microphone array to SH signals, based on a least-squares solution with
%   a constraint on noise amplification, using Tikhonov regularization. The
%   method formulates the LS problem in the spherical harmonic domain, by 
%   first expanding the array response to SH coefficients, and then solving
%   the LS problem, similar to
%
%       Jin, C.T., Epain, N. and Parthy, A., 2014. 
%       Design, optimization and evaluation of a dual-radius spherical microphone array. 
%       IEEE/ACM Transactions on Audio, Speech, and Language Processing, 22(1), pp.193-204.
%
%   Inputs:
%       H_array:    nFFT x nMics x nMeasurement_dirs matrix of measured
%           responses, in the frequency domain with length nFFT
%       order_sht:  order of SH signals to generate
%       grid_dirs:  nMeasurement_dirs x 2 matrix of [azi elev] of the
%           measurement directions, in rads
%       w_grid:     nMeasurement_dirs x 1 vector of weights for
%           weighted-least square solution, based on the importance or area
%           around each measurement point (leave empty if not known, or
%           not important)
%       nFFT:       number of FFT points in the responses
%       amp_threshold:  max allowed amplification for filters, in dB
%
%   Outputs:
%       H_filt: nSH x nMics x (nFFT/2+1) returned filters in the frequency 
%           domain (half spectrum up to Nyquist)
%       h_filt: nFFT x nMics x nSH impulse responses of the above filters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ARRAYSHTFILTERSMEAS_REGLSHD.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nGrid = length(H_array,3);
nMics = size(H_array,2);
if order_sht>sqrt(nMics)-1
    warning('Set order too high for the number of microphones, should be N<=sqrt(Q)-1')
    order_sht = floor( sqrt(nMics)-1 );
end
nBins = nFFT/2+1;
if isempty(w_grid)
    w_grid = ones(nGrid,1);
end

% SH matrix at grid directions
order_array = floor(sqrt(nGrid)/2-1);
aziElev2aziPolar = @(dirs) [dirs(:,1) pi/2-dirs(:,2)]; % function to convert from azimuth-inclination to azimuth-elevation
Y_grid = sqrt(4*pi) * getSH(order_array, aziElev2aziPolar(grid_dirs_rad), 'real')'; % SH matrix for grid directions

% compute inverse matrix
a_dB = amp_threshold;
alpha = 10^(a_dB/20);
beta = 1/(2*alpha);
W_grid = diag(w_grid);
H_filt = zeros((order_sht+1)^2, nMics, nBins);
% compute first the SHT of the array response
for kk=1:nBins
    tempH = squeeze(H_array(kk,:,:));
    H_nm(kk,:,:) = tempH * W_grid * Y_grid'* inv(Y_grid*W_grid*Y_grid');
end
% compute the inverse matrix in the SHD with regularization
for kk=1:nBins
    tempH_N = squeeze(H_nm(kk,:,:));
    tempH_N_trunc = tempH_N(:,1:(order_sht+1)^2);
    H_filt(:,:,kk) = tempH_N_trunc' * inv(tempH_N*tempH_N' + beta^2*eye(nMics));
end

if nargout>1
    % time domain filters
    h_filt = H_filt;
    h_filt(:,:,end) = abs(h_filt(:,:,end));
    h_filt = cat(3, h_filt, conj(h_filt(:,:,end-1:-1:2)));
    h_filt = real(ifft(h_filt, [], 3));
    h_filt = fftshift(h_filt, 3);
end
