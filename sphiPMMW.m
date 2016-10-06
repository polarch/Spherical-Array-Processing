function [W_pmmw, Pd_est, Ps_est] = sphiPMMW(Phi_x, Phi_n, src_dirs)
%SPHIPMMW Get beamforming weights for the iPMMW beamformer in the SHD
%   
%   Computes the beamforming weights of the informed Multi-channel 
%   Multi-Plane wave spatial filter (iPMMW), as introduced in 
%
%       Thiergart, O., Taseska, M. and Habets, E.A., 2014. 
%       An informed parametric spatial filter based on instantaneous direction-of-arrival estimates. 
%       IEEE/ACM Transactions on Audio, Speech, and Language Processing, 22(12), pp.2182-2196.   
%
%   but in the SHD instead of the microphone signal domain.
%
%   Inputs:
%       Phi_x:  (order+1)^2x(order+1)^2 covariance/correlation matrix
%           of SH signals
%       Phi_n:  (order+1)^2x(order+1)^2 covariance/correlation matrix
%           of noise in SH signals, estimated before
%       src_dirs:   Kx2 [azi elev] estimated or known DoAs of source signals 
%           in the sound scene 
%
%   Outputs:
%       W_pmmw: (order+1)^2xK matrix of beamforming weights, one set per
%           DoA
%       Pd_est: estimate of the diffuse sound power
%       Ps_est: estimate of the Kx1 source signal powers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHIPMMW.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSH = size(Phi_x,1);
order = sqrt(nSH)-1;
nSrc = size(src_dirs,1);

% steering vectors for source directions
src_dirs2 = [src_dirs(:,1) pi/2-src_dirs(:,2)]; % from azi-elevation to azi-inclination
Y_src = getSH(order, src_dirs2, 'real');
A = Y_src';

% estimate diffuse PSD
AA = A*A';
V_AA = sorted_eig(AA); % EVD on the steering vectors
U_AA = V_AA(:,nSrc+1:end); % eigenvectors of nullspace
D = U_AA'*Phi_n*U_AA;
E = U_AA'*U_AA;
[V_DE, S_DE] = eig(D,-E); % generalized EVD
[~,idx] = max(diag(S_DE));
c_d = V_DE(:,idx); % generalized eigenvector that maximizes diffuse-to-noise ratio
w_d = U_AA*c_d;
Pd_est = w_d'*Phi_x*w_d/(w_d'*w_d);

% estimate source signals PSD
for ns=1:nSrc
    C_ns = A(:,ns)*A(:,ns)';
    D_src(:,ns) = C_ns(:);
end
Phi_u = (Phi_n + Pd_est*eye(nSH)/(4*pi));
Phi_v = Phi_x - Phi_u;
Ps_est = (D_src'*D_src)\D_src' * Phi_v(:);
Phi_s = diag(Ps_est);

% trade-off beamforming control parameters
a_1 = 3;
a_2 = 0.5;
a_3 = -28;
sig = @(ksi) a_1/2 * ( 1 + tanh( a_2*(a_3-ksi)/2 ) );
phi_u = Pd_est + trace(Phi_n)/nSH;
ksi_l = 10*log10(Ps_est/phi_u);
beta_l = sig(ksi_l);
B = diag(beta_l);

% beamforming weights
W_pmmw = Phi_u\ A / (B/Phi_s + A'/Phi_u * A);

end
