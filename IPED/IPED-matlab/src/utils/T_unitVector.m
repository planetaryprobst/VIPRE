function [T_uv] = T_unitVector(S_uv)
%--------------------------------------------------------------------------
%   Returns the unit vector of the T-axis of the B-plane for an incoming
%   hyperbolic entry trajectory.
%--------------------------------------------------------------------------
%   Form:
%   [T_uv] = T_unitVector(S_uv)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   S_uv        (3,:)    -          unit vector S
%
%   ------
%   Output
%   ------
%   T_uv        (3,:)    -          unit vector T-axis of B-plane
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 12.01.2020 |  A. Probst    | First revision
%*************************************************************************%

% z-axis in planet's reference frame
n = [0;0;1];

% T_uv
T_uv = cross(S_uv,n)/vecnorm(cross(S_uv,n),2,1);


end



