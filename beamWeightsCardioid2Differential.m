function a_n = beamWeightsCardioid2Differential(N)
%BEAMWEIGHTSCARDIOID2DIFFERENTIAL Generate beamweights for differential array
%(cardioid)
%
%   For a specific order N of a higher order cardioid of the form
%   D(theta)=(1/2)^N * (1+cos(theta))^N, generate the beamweights for a
%   differential array of the form 
%   D(theta)=Sum_{n=1}^{N} a_n * cos^n(theta), that result in the cardioid.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSCARDIOID2DIFFERENTIAL.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The coefficients are given directly by the binomial expansion of the
% cardioid form
a_n = zeros(N+1,1);
for n=0:N
    a_n(n+1) = (1/2)^N * nchoosek(N, n);
end

end

