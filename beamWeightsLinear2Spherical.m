function b_n = beamWeightsLinear2Spherical(a_n, PLOT_ON)
%BEAMWEIGHTSLINEAR2SPHERICAL Convert linear array weights to spherical beamweights.
%
%   For a set of N+1 real a_n weights of the sensors of a 2*N+1 symmetrically 
%   weighted array, generate the beamweights for the same pattern in the SHD. 
%   Since the pattern is axisymmetric only the N+1 coefficients of m=0 are 
%   returned.Implemented from 
%
%       Hafizovic, I., Nilsen, C.I.C. and Holm, S., 2012. 
%       Transformation between uniform linear and spherical microphone arrays with symmetric responses. 
%       IEEE Transactions on Audio, Speech, and Language Processing, 20(4), pp.1189-1195.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSLINEAR2SPHERICAL.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% order N (or polynomial degree)
N = length(a_n)-1;
% linear pattern of symmetrically real-weighted 2*N+1 elements array, in terms of cosines
W_lin = eye(N+1); 
W_lin(2:end,2:end) = 2*W_lin(2:end,2:end); %weights for the coefficients

% build Chebyshev polynomial matrix
T_N = zeros(N+1);
for n = 0:N
    T_N(1:n+1,n+1) = chebyshevPoly(n);
end
% build Legendre polynomial matrix
P_N = zeros(N+1);
for n = 0:N
    P_N(1:n+1,n+1) = legendrePoly(n);
end

% build the normalisation matrix for the spherical array
w_sph = sqrt((2*(0:N)+1)/(4*pi));
W_sph = diag(w_sph);
W_sph_inv = diag(1./w_sph);

% spherical coefficients with respect to the linear array ones
b_n = W_sph_inv * inv(P_N) * T_N * W_lin * a_n;

% plot response if required
if PLOT_ON

    % cosine basis for linear array response
    theta = (0:1:360)'*pi/180;
    c_N = [];
    for n=0:N
        c_N = [c_N cos(n*pi*cos(theta))];
    end
    f_lin = c_N * W_lin * a_n;
    polar(theta, abs(f_lin))
    
    % build the legendre base
    p_N = [];
    theta2 = pi*cos(theta);    
    for n=0:N
        temp_leg = legendre(n, cos(theta2));
        p_N = [p_N temp_leg(1,:)'];
    end

    f_sph = p_N * W_sph * b_n;
    hold on, polar(theta, abs(f_sph), '--r')
    legend('linear array','spherical array')
end
