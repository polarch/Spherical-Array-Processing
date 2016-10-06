function diff = getDiffuseness_DPV(SHcov)
%GETDIFFUSENESS_DPV Compute the directional power variation diffuseness 
%estimator from SH signals
%   
%   This routine estimates diffuseness using the directional diffuseness of
%   the directional energy variance, as introduced by 
%
%       Thiele, R., 1953. Richtungsverteilung und zeitfolge der schallr?ckw?rfe in r?umen. 
%       Acta Acustica United with Acustica, 3(Supplement 2), pp.291-302.
%
%   and reformulated for spherical arrays by
%
%       Gover, B.N., Ryan, J.G. and Stinson, M.R., 2002. 
%       Microphone array measurement system for analysis of directional and spatial variations of sound fields. 
%       The Journal of the Acoustical Society of America, 112(5), pp.1980-1991.
%
%   This estimator is simplified somewhat here exploiting SHD properties. 
%   It is a robust estimator to the presence of multiple uncorrelated 
%   sources, but not to the presence of correlated ones, as it will still 
%   result in high diffuseness even if the sources are completely coherent.
%
%   Inputs:
%       SHcov:  (order+1)^2x(order+1)^2 covariance/correlation matrix of 
%       the SH signals
%
%   Outputs:
%       diff:       diffuseness value
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GETDIFFUSENESS_DPV.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aziElev2aziPolar = @(dirs) [dirs(:,1) pi/2-dirs(:,2)]; % function to convert from azimuth-inclination to azimuth-elevation

nSH = size(SHcov,1);
order = sqrt(nSH)-1;

% get grid for power estimation
[~,grid_dirs] = getTdesign(4*order);
Ygrid = getSH(order, aziElev2aziPolar(grid_dirs), 'real');
Agrid = (4*pi)/(order+1)^2*Ygrid;

% get directional power
p_grid = real(diag(Agrid*SHcov*Agrid'));
% mean directional power
p_mean = mean(p_grid);
% power deviation
p_dev = (1/p_mean)*sum( abs(p_grid - p_mean) );

% directivity factor for plane-wave decomposition (and max-DI) beamformer
pw_df = (order+1)^2;
y0 = getSH(order, [0 0], 'real')'; % unit power plane wave from 0,0
pw_grid = abs(Agrid*y0).^2;
pw_dev = pw_df*sum( abs(pw_grid - 1/pw_df) );

% diffuseness
diff = 1 - p_dev/pw_dev;

end
