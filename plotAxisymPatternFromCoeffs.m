function [ h_ax ] = plotAxisymPatternFromCoeffs(b_n, h_ax)
%PLOTAXISYMPATTERNFROMCOEFFS Plot axisymmetric pattern from SHD beamweights
%
%   For a specific order N, plot an axisymmetric beampattern described by 
%   the N+1 beamweights in the SHD, of m=0. The optional axis handle h_ax 
%   can be used to plot to a specified axis.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PLOTAXISYMPATTERNFROMCOEFFS.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    figure
    h_ax = gca;
end
    % build the SH basis for plotting (keep only the m=0 terms)
    theta = (0:1:360)'*pi/180;
    N = length(b_n)-1;
    y_nm = getSH(N, [zeros(size(theta)) theta], 'real');
    for n=0:N
        y_n(:,n+1) = y_nm(:,n^2+n+1);
    end
    % evaluate beampattern at specified directions
    f = y_n * b_n;
    
    axes(h_ax)
    polar(theta(f>=0), abs(f(f>=0)), 'b')
    hold on, polar(theta(f<0), abs(f(f<0)), 'r')
end

