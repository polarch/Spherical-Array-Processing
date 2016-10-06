function w_nulldiff = sphNullformer_diff(order, src_dirs)
%SPHNULLFORMER_DIFF Beamweights for a beamformer with nulls at specified directions, 
%and an omnidirectional constraint
%   
%   For a set of K directions, computes the beamforming weights that places
%   nulls at them while trying to keep the overall response as
%   omnidirectional as possible. Useful when trying to extract
%   reverberant/diffuse sound while suppressing strong directional sources.
%   From a statistical formulation, it is similar to the beamformer
%   proposed in
%
%       Thiergart, O. and Habets, E.A., 2014. 
%       Extracting reverberant sound using a linearly constrained minimum variance spatial filter. 
%       IEEE Signal Processing Letters, 21(5), pp.630-634.
%
%   but in the SHD instead of the microphone signal domain.
%
%   Inputs:
%       order:  order of SH signals
%       src_dirs:  Kx2 [azi elev] directions to place the nulls
%
%   Outputs:
%       w_nulldiff:  (order+1)^2xK matrix of beamweights, one set per
%           column for each beamforming direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHNULLFORMER_DIFF.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nSH = (order+1)^2;
    nSrc = size(src_dirs,1);
    src_dirs2 = [src_dirs(:,1) pi/2-src_dirs(:,2)]; % convert from azi-elev to azi-incl

    % steering vectors in the SHD
    Y_nulldiff = getSH(order, src_dirs2, 'real');

    % compute weight vector
    c = zeros(nSrc+1,1);
    c(1) = 1; % constraint vector    
    g = zeros(nSH,1);
    g(1) = 1;
    A = [g Y_nulldiff']; 
    w_nulldiff = pinv(A')*c;
end
