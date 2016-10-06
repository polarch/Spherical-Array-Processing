function T_n = returnChebyPolyCoeffs(n)
%RETURNCHEBYPOLYCOEFFS Generate polynomial coefficients of a Chebyshev polynomial.
%
%   Generate the polynomial coefficients of a Chebyshev polynomial of 
%   degree n, and return them as a vector.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RETURNCHEBYPOLYCOEFFS.M - 10/4/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n==0
    T_n = 1;
else
    
    for k=0:floor(n/2)
        coeffs(k+1) = (n/2)*(-1)^k * factorial(n-k-1)/(factorial(k)*factorial(n-2*k)) * 2^(n-2*k);
    end
    
    T_n = zeros(n+1,1);
    if mod(n,2)==0
        T_n(1:2:end) = coeffs(end:-1:1);
    else
        T_n(2:2:end) = coeffs(end:-1:1);
    end
end
