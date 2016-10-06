function diff = getDiffuseness_IE(pvCOV)
%GETDIFFUSENESS_IE Compute the intensity-energy density ratio diffuseness 
% estimator from pressure-velocity (omni-dipole) signals
%   
%   This routine estimates diffuseness using the acoustic intensity - energy 
%   density ratio, as popularized by the Spatial Impulse Response Rendering 
%   (SIRR) and Directional Audio Coding (DirAC) methods. The estimator is
%   based on the relation 1 - ||Re{I_a}||/(c*E), where Ia is the active
%   intensity vector, and E the energy density. All the quantities can be 
%   computed from the pressure and velocity (dipole) signals, which in turn
%   can be computed from first-order SH signals.
%   This estimator gives correct results in the mixture of a point/plane 
%   wave source and isotropic diffuse sound, but over-estimates
%   significantly diffuseness for two incoherent or coherent sources, also
%   dependent on their directional separation.
%
%   Inputs:
%       PVcov:  4x4 covariance/correlation matrix of pressure (omni) and
%       three velocity (dipole) signals, in [p v_x v_y v_z] order
%
%   Outputs:
%       diff:       diffuseness value
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GETDIFFUSENESS_IE.M - 5/10/2016
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% active intensity vector
Ia = real(pvCOV(2:4,1));
Ia_norm = sqrt(sum(Ia.^2));
% energy estimate
E = trace(pvCOV)/2;

diff = 1 - Ia_norm/E;

end
