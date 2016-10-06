function P_n = returnLegePolyCoeffs(n)
%RETURNLEGEPOLYCOEFFS Generate polynomial coefficients of a Legendre polynomial.
%
%   Generate the polynomial coefficients of a Legendre polynomial of degree
%   n, and return them as a vector.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RETURNLEGEPOLYCOEFFS.M - 10/4/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=0:floor(n/2)
    coeffs(k+1) = (1/2^n)*(-1)^k * factorial(n)/(factorial(k)*factorial(n-k)) * factorial(2*n-2*k)/(factorial(n)*factorial(n-2*k));
end

P_n = zeros(n+1,1);
if mod(n,2)==0
    P_n(1:2:end) = coeffs(end:-1:1);
else
    P_n(2:2:end) = coeffs(end:-1:1);
end
