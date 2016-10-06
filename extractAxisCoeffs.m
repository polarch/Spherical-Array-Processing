function a_n = extractAxisCoeffs(a_nm)
%EXTRACTAXISCOEFFS Exctract the N+1 coeffs of m=0 from the full (order+1)^2
% coefficient vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXTRACTAXISCOEFFS.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = sqrt(size(a_nm,1))-1;
nSets = size(a_nm,2);
a_n = zeros(N+1,nSets);
for ns=1:nSets
    for n=0:N
        a_n(n+1,ns) = a_nm(n^2 + n + 1);
    end
end

end
