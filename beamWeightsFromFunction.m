function b_nm = beamWeightsFromFunction(fHandleArray, order)
%BEAMWEIGHTSFROMFUNCTION Generate spherical beamweights for arbitrary pattern
%
%   This function computes the SH coefficients of an arbitrary target
%   beampattern, by performing the SHT on it. The pattern should be
%   specified as a function handle that expects a Kx1 vector of azimuths as
%   first argument, and a vector of Kx1 elevations as the second argument, 
%   and return the functions values for these K directions.
%   Multiple patterns can be passed as an array of function handles.
%   The order specifies the maximum SH order of the beamweights.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSFROMFUNCTION.M - 10/4/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get a dense uniform design
[~, dirsAziElev] = getTdesign(20);

Nfunc = length(fHandleArray);
b_nm = zeros((order+1)^2, Nfunc);
for nf =1:Nfunc
    if Nfunc==1, fPattern = fHandleArray;
    else fPattern = fHandleArray{nf};
    end

    % compute the pattern at grid points
    F_dirs = fPattern(dirsAziElev(:,1), dirsAziElev(:,2));

    % get harmonic coefficients
    b_nm(:,nf) = directSHT(order, F_dirs, dirsAziInc, 'real');
end
