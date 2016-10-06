function [H_filt, h_filt] = arraySHTfiltersTheory_radInverse(R, Nmic, order_sht, Lfilt, fs, amp_threshold)
%ARRAYSHTFILTERSTHEORY_RADINVERSE Generate SHT filters based on theoretical
%responses (regularized radial inversion)
%
%   Generate the filters to convert microphone signals from a spherical
%   microphone array to SH signals, based on an ideal theoretical model of 
%   the array. The filters are generated from an inversion of the radial 
%   components of the response, neglecting spatial aliasing effects and
%   non-ideal arrangements of the microphones. One filter is shared by all
%   SH signals of the same order in this case.
%   Here this single channel inversion problem is done with a constraint on
%   max allowed amplification regularization, using Tikhonov
%   regularization, similar e.g. to 
%
%       Moreau, S., Daniel, J., Bertet, S., 2006, 
%       3D sound field recording with higher order ambisonics-objective 
%       measurements and validation of spherical microphone. 
%       In Audio Engineering Society Convention 120.
%
%   Inputs:
%       R:          radius of spherical array
%       Nmic:       number of microphones on the sphere
%       order_sht:  order of SH signals to generate
%       Lfilt:      number of FFT points for the output filters
%       fs:         sample rate for the output filters
%       amp_threshold:      max allowed amplification for filters, in dB
%
%   Outputs:
%       H_filt: nSH x nMics x (Lfilt/2+1) returned filters in the frequency 
%           domain (half spectrum up to Nyquist)
%       h_filt: Lfilt x order_sht+1 impulse responses of the above filters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ARRAYSHTFILTERSTHEORY_RADINVERSE.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 343;
f = (0:Lfilt/2)'*fs/Lfilt;

if order_sht>sqrt(Nmic)-1
    warning('Set order too high for the number of microphones, should be N<=sqrt(Q)-1')
    order_sht = floor( sqrt(Nmic)-1 );
end

% modal responses
kR = 2*pi*f*R/c;
bN = sphModalCoeffs(order_sht, kR, 'rigid')/(4*pi);

% encoding matrix with regularization
a_dB = amp_threshold;
alpha = sqrt(Nmic)*10^(a_dB/20);
beta = sqrt((1-sqrt(1-1/alpha^2))/(1+sqrt(1-1/alpha^2))); % Moreau & Daniel
% regularized single channel equalization filters per order
H_filt = conj(bN)./(abs(bN).^2 + beta^2*ones(Lfilt/2+1,order_sht+1));

if nargout>1
    % time domain filters
    h_filt = H_filt;
    h_filt(end,:) = real(h_filt(end,:));
    h_filt = [h_filt; conj(h_filt(end-1:-1:2,:))];
    h_filt = real(ifft(h_filt));
    h_filt = fftshift(h_filt, 1);    
end