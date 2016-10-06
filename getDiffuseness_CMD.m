function [diff, diff_ord] = getDiffuseness_CMD(SHcov)
%GETDIFFUSENESS_CMD Compute the COMEDIE diffuseness estimator from SH signals
%   
%   This routine estimates diffuseness by the COMEDIE estimator as
%   introduced in 
%
%       Epain, N. and Jin, C.T., 2016. 
%       Spherical Harmonic Signal Covariance and Sound Field Diffuseness. 
%       IEEE/ACM Transactions on Audio, Speech, and Language Processing, 24(10), pp.1796-1807.
%
%   This estimator is the most advanced from the ones included in the
%   library, as it takes into account correlations between higher-order 
%   signals if available, rises gradually in the presence of multiple 
%   uncorrelated sources, and is robust in the presence of correlated ones.
%
%   Inputs:
%       SHcov:  (order+1)^2x(order+1)^2 covariance/correlation matrix of 
%       the SH signals
%
%   Outputs:
%       diff:       diffuseness value
%       diff_ord:   diffuseness value for each of the lower orders in the
%           signals, if asked
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GETDIFFUSENESS_CMD.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSH = size(SHcov,1);
order = sqrt(nSH)-1;

%eigendecomposition - sort eigenvalues in descending order
[~, Dsh] = eig(SHcov);

%diffuseness estimation
g_0 = 2*((order+1)^2-1);
mean_ev = (1/(order+1)^2)*sum(diag(Dsh));
g = (1/mean_ev)*sum(abs(diag(Dsh)-mean_ev));
diff = 1-g/g_0;

%if it is called, return diffuseness per order
if nargout>1
    diff_ord = zeros(1,order);
    for n=1:order-1
        %covariance matrices
        SHcov_n = SHcov(1:(n+1)^2,1:(n+1)^2);
        
        %eigendecomposition - sort eigenvalues in descending order
        [~, Dsh_n] = eig(SHcov_n);
        
        %diffuseness estimation
        g_n0 = 2*((n+1)^2-1);
        mean_ev = (1/(n+1)^2)*sum(diag(Dsh_n));
        g_n = (1/mean_ev)*sum(abs(diag(Dsh_n)-mean_ev));
        diff_ord(n) = 1-g_n/g_n0;
        
    end
    diff_ord(order) = diff;
    
end

end
