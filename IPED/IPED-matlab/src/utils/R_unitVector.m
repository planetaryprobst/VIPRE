function [R_uv] = R_unitVector(S_uv,T_uv)
%--------------------------------------------------------------------------
%   Returns the unit vector of the R-axis 'R_uv' of the B-plane for an
%   incoming hyperbolic entry trajectory using the unit vector 'S_uv' and
%   'T_uv'.
%--------------------------------------------------------------------------
%   Form:
%   [R_uv] = R_unitVector(S_uv,T_uv)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   S_uv        (3,:)    -          unit vector S
%   T_uv        (3,:)    -          unit vector T-axis of B-plane
%
%   ------
%   Output
%   ------
%   R_uv        (3,:)    -          unit vector R-axis of B-plane
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 12.01.2020 |  A. Probst    | First revision
%*************************************************************************%

% T_uv
R_uv = cross(S_uv,T_uv);

end



