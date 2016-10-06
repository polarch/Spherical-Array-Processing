function d_n = beamWeightsDolphChebyshev2Spherical(N, paramType, arrayParam )
%BEAMWEIGHTSDOLPHCHEBYSHEV2SPHERICAL Generate spherical coefficients for
%Dolph-Chebyshev designs
%
%   Generate beamweights in the SHD for Dolph-Chebyshev beampatterns, with
%   mainlobe and sidelobe control. Because the pattern is axisymmetric only 
%   the N+1 coefficients of m=0 are returned.
%
%   N:          order
%   paramType:  either 'sidelobe', for sidelobe leve control, or 'width'
%       for mainlobe width control
%   arrayParam: either sidelobe level 1/R, or mainlobe width 2*a0
%
%   Based on the implementation of
%
%       Koretz, A. and Rafaely, B., 2009. 
%       Dolph-Chebyshev beampattern design for spherical arrays. 
%       IEEE Transactions on Signal Processing, 57(6), pp.2417-2420.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSDOLPHCHEBYSHEV2SPHERICAL.M - 10/4/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% order of Chebyshev polynomial (pattern)
M = 2*N;

if isequal(paramType, 'sidelobe')
    R = 1/arrayParam;
    x0 = cosh((1/M) * acosh(R));
elseif isequal(paramType, 'width')
    a0 = arrayParam/2;
    x0 = cos(pi/(2*M))/cos(a0/2);
    R = cosh(M*acosh(x0));
end

% compute Chebyshev polynomial vector of order 2N
t_2N = returnChebyPolyCoeffs(2*N);

% build Legendre polynomial matrix up to order N
P_N = zeros(N+1);
for n = 0:N
    P_N(1:n+1,n+1) = returnLegePolyCoeffs(n);
end

% compute coefficients of the Dolph-Chebyshev pattern
d_n=zeros(N+1,1);
for n=0:N

    temp = 0;
    for i=0:n
        for j=0:N
            for m=0:j
                temp = temp+(1-(-1)^(m+i+1))/(m+i+1)*factorial(j)/(factorial(m)*factorial(j-m))*(1/2^j)*t_2N(2*j+1)*P_N(i+1,n+1)*x0^(2*j);
            end
        end
    end
    d_n(n+1) = (2*pi/R)*temp;    
end

% normalize coefficients to unity response at look-direction
norm = sum(d_n .* sqrt(2*(0:N)+1)' );
d_n = sqrt(4*pi) * d_n/norm;
