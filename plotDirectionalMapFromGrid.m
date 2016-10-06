function h_ax = plotDirectionalMapFromGrid(fgrid, aziRes, polarRes, h_ax, POLAR_OR_ELEV, ZEROED_OR_CENTERED)
%PLOTDIRECTIONALMAPFROMGRID Compute the intensity-energy density ratio diffuseness 
% estimator from pressure-velocity (omni-dipole) signals
%   
%   Plots a directional function evluated on a regular grid of
%   azimuth-elevation intervals. The function values are passed as a vector
%   with ordering according to what the function grid2dirs.m (SHT-lib) returns.
%
%   Inputs:
%       fgrid:  Kx1 vector of function values to be plotted on the grid
%       aziRes: resolution in degrees for azimuth, should divide 360
%           exactly and it should be the same as used in grid2dirs.m
%       polarRes:   resolution in elevation in degrees, should divide 180
%           exactly and it should be the same as used in grid2dirs.m
%       h_ax:   optional axis handle to use, otherwise create new
%       POLAR_OR_ELEV:  flag for directional indexing, 1 if the polar
%           angles denote inclinations from the top (0~180), or 0 if they
%           denote elevations from the horizon (-90~90, Matlab
%           convention). It should be the same as used in grid2dirs.m.
%           Default is 1.
%       ZEROED_OR_:  flag for directional indexing, 1 if the azimuth angles
%           range from 0~360, or 0 if they range from -180~180. It should 
%           be the same as used in grid2dirs.m. Default is 1.
%
%   Outputs:
%       h_ax:   axis handle of the axis that was used or generated
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PLOTDIRECTIONALMAPFROMGRID.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 5
        ZEROED_OR_CENTERED = 1;
    case 4
        POLAR_OR_ELEV = 1;
        ZEROED_OR_CENTERED = 1;
    case 3
        POLAR_OR_ELEV = 1;
        ZEROED_OR_CENTERED = 1;
        figure
        h_ax = gca;
    case {2,1,0}
        error('Not enough arguments')
end
if isempty(h_ax) || ~ishandle(h_ax)
    figure;
    h_ax = gca;
end
axes(h_ax);

Fgrid = Fdirs2grid(fgrid, aziRes, polarRes, 1);

if POLAR_OR_ELEV
    elev = 0:polarRes:180;
else
    elev = -90:polarRes:90;
end

if ZEROED_OR_CENTERED
    azi = 0:aziRes:360;
else
    azi = -180:aziRes:180;
end

if POLAR_OR_ELEV
    imagesc(azi,elev, Fgrid)
else
    imagesc(azi,elev, Fgrid), set(gca,'YDir','normal')
end

