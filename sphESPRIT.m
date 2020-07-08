function src_dirs_rad = sphESPRIT2(Us)
%SPHESPRIT DoA estimation using ESPRIT in the SHD
%   
%   This routine spherical ESPRIT method to SH signals and returns the 
%   analyzed DoAs without a grid search. It implements the 3-recurrence
%   relationship variant proposed by Jo & Choi in
%
%   B. Jo and J.-W. Choi, “Parametric direction-of-arrival estimation with
%   three recurrence relations of spherical harmonics,” 
%   J. Acoust. Soc. Amer.,vol. 145, no. 1, pp. 480–488, Jan. 2019.
%
%   Subspace approaches such as ESPRIT can offer higher spatial resolution 
%   than beamforming approaches such as the steered-response power, as long
%   as the source signals are not correlated between them and with the 
%   reverberant/diffuse sound.
%
%   Inputs:
%       Us: (order+1)^2xK signal subspace taken from the first K eigenvectors 
%       of the spatial correlation matrix, after sorting in eigenvalue 
%       descending order
%
%   Outputs:
%       est_dirs:   nSrcx2 [azi elev] of estimated directions from
%           peak-finding
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHESPRIT.M - 5/3/2019
% Archontis Politis, archontis.politis@tuni.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[LambdaXYp, LambdaXYm, LambdaZ] = getLambda(Us);
[PsiXYp, PsiXYm, PsiZ] = getPsi(Us, LambdaXYp, LambdaXYm, LambdaZ);
[V, ~] = eig(PsiXYp, PsiZ, 'qz'); % A*V = B*V*D -> [V,D]=eig(A,B).
PhiXYp = V\(PsiXYp*V);
PhiXYm = V\(PsiXYm*V);
PhiZ   = V\(PsiZ*V);
phiX   = real(diag(PhiXYp+PhiXYm)/2);
phiY   = real(diag(PhiXYp-PhiXYm)/(2i));
phiZ   = real(diag(PhiZ));

azim = atan2(phiY,phiX);
elev = atan2(phiZ,sqrt(phiX.^2+phiY.^2));
src_dirs_rad = [azim elev];

end


function Ynimu = getYnimu(Ynm, ni, mu)

N = sqrt(size(Ynm,2))-1;
[idx_nimu, idx_nm] = muni2q(N,ni,mu);
Ynimu = zeros(size(Ynm,1),N^2);
Ynimu(:,idx_nimu) = Ynm(:,idx_nm);

end


function [idx_nimu, idx_nm] = muni2q(order,ni,mu)

nm = [];
for n=0:order-1
    nm = [nm; n*ones(2*n+1,1) (-n:n)'];
end
nimu = [nm(:,1)+ni nm(:,2)+mu];
qnm = nm(:,1).^2+nm(:,1)+nm(:,2)+1;
qnimu = nimu(:,1).^2+nimu(:,1)+nimu(:,2)+1;
idx_valid = find(abs(nimu(:,2))<=nimu(:,1));
idx_nm = qnimu(idx_valid);
idx_nimu = qnm(idx_valid);

end


function Wnimu = getWnimu(order, mm, ni, mu)

nm = [];
for n=0:order-1
    nm = [nm; n*ones(2*n+1,1) (-n:n)'];
end
if mm==1
    nimu = [nm(:,1)+ni nm(:,2)+mu];
elseif mm==-1
    nimu = [nm(:,1)+ni -nm(:,2)+mu];
end
w_nimu = sqrt( (nimu(:,1)-nimu(:,2)-1).*(nimu(:,1)-nimu(:,2))./((2*nimu(:,1)-1).*(2*nimu(:,1)+1)) );
Wnimu = diag(w_nimu);

end


function Vnimu = getVnimu(order, ni, mu)

nm = [];
for n=0:order-1
    nm = [nm; n*ones(2*n+1,1) (-n:n)'];
end
nimu = [nm(:,1)+ni nm(:,2)+mu];
v_nimu = sqrt( (nimu(:,1)-nimu(:,2)).*(nimu(:,1)+nimu(:,2)) ./((2*nimu(:,1)-1).*(2*nimu(:,1)+1)) );
Vnimu = diag(v_nimu);

end


function [PsiXYp, PsiXYm, PsiZ] = getPsi(Us, LambdaXYp, LambdaXYm, LambdaZ)

pinvUs = pinv(getYnimu(Us.',0,0).');
PsiXYp = pinvUs*LambdaXYp;
PsiXYm = pinvUs*LambdaXYm;
PsiZ   = pinvUs*LambdaZ;

end


function [PhiXYp, PhiXYm, PhiZ] = getPhi(src_dirs_rad)

PhiZ = diag(cos(src_dirs_rad(:,2)));
PhiXYp = diag( sin(src_dirs_rad(:,2)).*exp(1i*src_dirs_rad(:,1)) );
PhiXYm = diag( sin(src_dirs_rad(:,2)).*exp(-1i*src_dirs_rad(:,1)) );

end


function [LambdaXYp, LambdaXYm, LambdaZ] = getLambda(Us)

order = sqrt(size(Us,1))-1;
LambdaXYp =  getWnimu(order, 1,1,-1)*getYnimu(Us.', 1,-1).' - getWnimu(order,-1,0,0)*getYnimu(Us.',-1,-1).';
LambdaXYm = -getWnimu(order,-1,1,-1)*getYnimu(Us.', 1, 1).' + getWnimu(order, 1,0,0)*getYnimu(Us.',-1, 1).';
LambdaZ   =  getVnimu(order,   0, 0)*getYnimu(Us.',-1, 0).' + getVnimu(order,   1,0)*getYnimu(Us.', 1, 0).';

end
