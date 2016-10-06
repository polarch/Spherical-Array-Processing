function gammaAB = diffCoherence(k, r_A, r_B, a_nm, b_nm, G_mtx)
%DIFFCOHERENCE Computes diffuse field coherence of two directional patterns
%
%   Does that. For more details see
%
%       Politis, A., 2016. 
%       Diffuse-field coherence of sensors with arbitrary directional responses. arXiv preprint arXiv:1608.07713.
%
%   Inputs:
%       k:  vector of wavenumbers to compute coherence at
%       r_A:    [x y z] vector of location of first sensor/beamformer
%       r_B:    [x y z] vector of location of second sensor/beamformer
%       a_nm:   SHD beamweights of the first directional pattern
%       b_nm:   SHD beamweights of the second directional pattern
%       G_mtx:  (optional) pre-computed matrix of Gaunt coefficients used
%           in the derivation. If not passed then they are computed here,
%           but that can be very slow for high orders.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DIFFCOHERENCE.M - 7/06/2015
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the vector connecting the two patterns
r_AB = r_B-r_A;
[phi_r, elev_r, d] = cart2sph(r_AB(1),r_AB(2),r_AB(3));
theta_r = pi/2-elev_r;

orderA = sqrt(length(a_nm)) -1;
orderB = sqrt(length(b_nm)) -1;
orderC = orderA+orderB;
minOrder = min(orderA, orderB);

% convert second pattern to its conjugate for coherence between complex
% patterns (doesn't do anything if the patterns are real)
b_nm = conjCoeffs(b_nm);

% denominator of the coherence function
denAB = sqrt((a_nm'*a_nm) * (b_nm'*b_nm));

if d==0
    % numerator of the coherence function
    numAB = dot(a_nm(1:(minOrder+1)^2), b_nm(1:(minOrder+1)^2));
    numAB = numAB*ones(size(k));
    
else
    % get coefficients of the product function
    if nargin<6
        c_nm = sphMultiplication(a_nm, b_nm);
    else
        c_nm = sphMultiplication(a_nm, b_nm, G_mtx);
    end
    
    % numerator of the coherence function
    y_nm_r = getSH(orderC, [phi_r theta_r], 'complex')';
    kd = k*d;
    imag_n = 1i.^(0:orderC).';
    bessel_n = zeros(orderC+1,length(kd));
    for no=0:orderC, bessel_n(no+1,:) = sph_besselj(no, kd).'; end
    imag_bessel_n = (imag_n*ones(1,length(kd))).*bessel_n;
    imag_bessel_nm = replicatePerOrder(imag_bessel_n);
    planewave_nm = 4*pi*(y_nm_r*ones(1,length(kd))).*imag_bessel_nm;
    numAB = zeros(size(k));
    for kk=1:length(kd)
        numAB(kk) = planewave_nm(:,kk)' * c_nm;
    end
end
    
% coherence
gammaAB = numAB/denAB;

end
