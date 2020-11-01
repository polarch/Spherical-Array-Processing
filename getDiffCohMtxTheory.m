function M_diffcoh = getDiffCohMtxTheory(mic_dirs_rad, arrayType, R, N_max, freqs, dirCoeff)
%GETDIFFCOHMTXTHEORY Returns the diffuse coherence matrix of atheoretical SMA
%   GETDIFFCOHMTXTHEORY computes the spatial covariance matrix of a 
%   theoretical spherical microphone array, by modeling its array response 
%   directly in the spherical harmonic domain, for a theoretical isotropic 
%   diffuse field of unit power. We refer to this as the diffuse field coherence 
%   matrix D of the array where each element (i,j) of the matrix represents 
%   the complex coherence between channel i and j, for diffuse excitation. 
%   Technically, the element should be also normalized with the signal powers 
%   to correspond to the standard coherence formula, as: D_ij/sqrt(DII*Djj). 
%   This matrix is very useful for diffuse field modeling and synthesis, 
%   diffuse-field equalization techniques, and spatial filtering that involves 
%   models of ambient noise as diffuse, for suppresion.
%
%   Inputs:
%       mic_dirs_rad :  [mic_azi1 mic_elev1; ...] directions of microphone
%           capsules
%       arrayType    :  'open', 'rigid', or 'directional'
%       R            :  array radius in m
%       N_max        :  maximum order of spherical approximation
%       fs           :  sample rate
%       dirCoeff     :  if array consists of directional microphones, then
%           dirCoeff is the directivity coefficient from 0 to 1 
%           (1 for omni, 0.5 cardioid, 0 dipole)
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

mic_dirs_xyz = unitSph2cart(mic_dirs_rad);
nMics = size(mic_dirs_xyz,1);
nFreqs = length(freqs);

c = 343;
bn = sphModalCoeffs(N_max, 2*pi*freqs*R/c, arrayType, dirCoeff)/(4*pi); % if 4pi is inside bn
bn2 = abs(bn).^2;

M_diffcoh = zeros(nMics,nMics,nFreqs);
for col=1:nMics
    for row=col:nMics
        cosangle = sum(mic_dirs_xyz(row,:).*mic_dirs_xyz(col,:));
        if cosangle>1, cosangle = 1; end
        if cosangle<-1, cosangle = -1; end
        for order=0:N_max
            temp = legendre(order,cosangle);
            Pn(order+1) = temp(1);
        end
        M_diffcoh(row,col,:) = 4*pi* bn2 * ((2*(0:N_max)+1) .* Pn).';
        M_diffcoh(col,row,:) = M_diffcoh(row,col,:);
    end
end

end
