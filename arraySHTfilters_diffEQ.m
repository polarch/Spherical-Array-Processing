function M_mic2sh_diffeq = arraySHTfilters_diffEQ(M_mic2sh, M_dfc, f_alias, fs)
%ARRAYSHTFILTERS_DIFFEQ Diffuse-filed equalization of the SMA encoder
%   ARRAYSHTFILTERS_DIFFEQ computes the diffuse-field energy of the aliased 
%   encoded components, and equalizes them to have flat diffuse-field response.
%   Note that this operation does not fix aliasing itself, it is just a
%   simple operation that corrects the energy of the aliased signals "on
%   average". It is also the approach above aliasing that was proposed by
%   Gerzon in the original Soundfield microphone publications.
%
%   Inputs:
%       M_mic2sh:   nSH x nMics x nBins matrix of spherical harmonic
%           encoding filters
%       M_dfc:      nMics x nMics x nBins diffuse-field coherence matrix of
%           the spherical microphone array (SMA)
%       f_alias:    some estimated aliasing frequency limit for the SMA
%       fs:         sampling frequency
%
%   Outputs:
%       M_mic2sh_diffeq:    nSH x nMics x nBins matrix of spherical harmonic
%           encoding filters, diffuse-field equalized above aliasing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ARRAYSHTFILTERS_DIFFEQ.M - 8/6/2017
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nBins   = size(M_mic2sh,3);
f = (0:nBins-1)*fs/(2*(nBins-1));

[~,idxf_alias] = min(abs(f-f_alias(1))); % find aliasing frequency
M_mic2sh_diffeq(:,:,1:idxf_alias) = M_mic2sh(:,:,1:idxf_alias);
E_fal = M_mic2sh(:,:,idxf_alias);
L_diff_fal = diag(E_fal*M_dfc(:,:,idxf_alias)*E_fal'/(4*pi));
for nf=idxf_alias+1:length(f)
    E = M_mic2sh(:,:,nf);
    L_diff_f = diag(E*M_dfc(:,:,nf)*E'/(4*pi));
    M_mic2sh_diffeq(:,:,nf) = diag(sqrt(L_diff_fal./L_diff_f)) * E;
end

end
