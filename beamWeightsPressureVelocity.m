function M_sh2pv = beamWeightsPressureVelocity(basisType)
% BEAMWEIGHTSPRESSUREVELOCITY Convert from SH signals to pressure velocity
%
%   Returns the 4x4 matrix that converts first-order SH signals to pressure
%   velocity signals, with the ordering [p v_x v_y v_z], where v_xyz are the
%   cartesian components of the acoustic veclocity vector.
%
%   basisType:  'complex' or 'real' for the respective type of spherical
%       harmonics under use.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BEAMWEIGHTSPRESSUREVELOCITY.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch basisType
    case 'real'
        M_sh2pv = sqrt(4*pi)*[  1   0           0           0;
                                0   0           0           1/sqrt(3);
                                0   1/sqrt(3)   0           0;
                                0   0           1/sqrt(3)   0];
    case 'complex'
        M_sh2pv = sqrt(4*pi)*[  1   0           0           0;
                                0   1/sqrt(6)   0           -1/sqrt(6);
                                0   -1i/sqrt(6) 0           -1i/sqrt(6);
                                0   0           1/sqrt(3)   0];
end
