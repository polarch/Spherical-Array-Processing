function A_xyz = computeVelCoeffsMtx(sectorOrder)
%COMPUTEVELCOEFFSMTX Compute the generating matrices of products of
%beampatterns with dipoles
%
%   This routine computes the matrices that generate the coefficients of
%   the beampattern of order (sectorOrder+1) that is essentially the 
%   product of a pattern of order=sectorOrder and a dipole. It is used in 
%   beamWeightsVelocityPatterns.m. For the derivation of the matrices see:
%
%       Politis, A. and Pulkki, V., 2016. 
%       Acoustic intensity, energy-density and diffuseness estimation in a directionally-constrained region. 
%       arXiv preprint arXiv:1609.03409.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% COMPUTEVELCOEFFSMTX.M - 7/06/2015
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% orders of sector and resulting velocity patterns
Ns = sectorOrder;
Nxyz = Ns+1;

% Complex SH coefficients of velocity x,y,z
% X = cos\phi sin\theta
% Y = sin\phi sin\theta
% Z = cos\theta
% q = 1 -> Yq = Y_{1(-1)}
% q = 2 -> Yq = Y_{10}
% q = 3 -> Yq = Y_{11}
% x_q != 0 only for q=1,3
% y_q != 0 only for q=1,3
% z_q != 0 only for q=2
x1 = sqrt(2*pi/3);
x3 = -sqrt(2*pi/3);
y1 = 1i*sqrt(2*pi/3);
y3 = 1i*sqrt(2*pi/3);
z2  = sqrt(4*pi/3);

A_xyz = zeros((Nxyz+1)^2, (Ns+1)^2, 3);

G_mtx = gaunt_mtx(Ns, 1, Nxyz);

for ii=1:(Nxyz+1)^2
    for jj=1:(Ns+1)^2
        
        A_xyz(ii,jj,1) = x1*G_mtx(jj,2,ii) + x3*G_mtx(jj,4,ii);
        A_xyz(ii,jj,2) = y1*G_mtx(jj,2,ii) + y3*G_mtx(jj,4,ii);
        A_xyz(ii,jj,3) = z2*G_mtx(jj,3,ii);
    end
end
