
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSDIFFERENTIAL2SPHERICAL.M - 10/4/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tabulated coefficients of the differential form for common beampatterns
% up to fourth-order, copied from
%
%       Elko, G.W., 2004. Differential microphone arrays. 
%       In Audio signal processing for next-generation multimedia communication systems (pp. 11-65). Springer.
%

% Cardioids
a1 = [1/2 1/2]';
a2 = [1/4 2/4 1/4]';
a3 = [1/8 3/8 3/8 1/8]';
a4 = [1/16 4/16 6/16 4/16 1/16]';

% Supercardioids
a1 = [sqrt(3)-1 3-sqrt(3)]'/2;
a2 = [1 2*sqrt(7) 5]'/2/(3+sqrt(7));

b = sqrt(2);
c = sqrt(21);
d = sqrt(21-c);
a3 = [b*d-c-1  21+9*c-b*(6+c)*d   3*(b*(4+c)*d-25-5*c)  63+7*c-b*(7+2*c)*d]'/8;
a4 = [0.0036 0.0670 0.2870 0.4318 0.2107]';

% Hypercardioids
a1 = [1/4 3/4]';
a2 = [-1/6 2/6 5/6]';
a3 = [-3/32 -15/32 15/32 35/32]';
a4 = [0.075 -0.3 -1.05 0.7 1.575]';
 