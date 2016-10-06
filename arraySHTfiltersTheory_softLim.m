function [H_filt, h_filt] = arraySHTfiltersTheory_softLim(R, Nmic, order_sht, Lfilt, fs, amp_threshold)
%ARRAYSHTFILTERSTHEORY_SOFTLIM Generate SHT filters based on theoretical
%responses (soft-limiting)
%
%   Generate the filters to convert microphone signals from a spherical
%   microphone array to SH signals, based on an ideal theoretical model of 
%   the array. The filters are generated from an inversion of the radial 
%   components of the response, neglecting spatial aliasing effects and
%   non-ideal arrangements of the microphones. One filter is shared by all
%   SH signals of the same order in this case.
%   Here this single channel inversion problem is done through a
%   thresholding approach on the inversion, limited to a max allowed
%   amplification. The limiting follows the approach of 
%
%       Bernschutz, B., Porschmann, C., Spors, S., Weinzierl, S., Versterkung, B., 2011. 
%       Soft-limiting der modalen amplitudenverst?rkung bei sph?rischen mikrofonarrays im plane wave decomposition verfahren. 
%       Proceedings of the 37. Deutsche Jahrestagung fur Akustik (DAGA 2011)
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
%       H_filt: (Lfilt/2+1) x order_sht+1 returned filters in the frequency 
%           domain (half spectrum up to Nyquist)
%       h_filt: Lfilt x order_sht+1 impulse responses of the above filters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ARRAYSHTFILTERSTHEORY_SOFTLIM.M - 5/10/2016
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
bN = sphModalCoeffs(order_sht, kR, 'rigid')/(4*pi); % due to modified SHs, the 4pi term disappears from the plane wave expansion
% single channel equalization filters per order
inv_bN = 1./bN;
inv_bN(1,2:end) = 0; % omit INFs at DC

% encoding matrix with soft-limiting
a_dB = amp_threshold;
alpha = sqrt(Nmic)*10^(a_dB/20);
% regularized single channel equalization filters per order
H_filt = (2*alpha/pi) * (abs(bN).*inv_bN).*atan((pi/(2*alpha))*abs(inv_bN));

if nargout>1
    % time domain filters
    h_filt = H_filt;
    h_filt(end,:) = real(h_filt(end,:));
    h_filt = [h_filt; conj(h_filt(end-1:-1:2,:))];
    h_filt = real(ifft(h_filt));
    h_filt = fftshift(h_filt, 1);
end