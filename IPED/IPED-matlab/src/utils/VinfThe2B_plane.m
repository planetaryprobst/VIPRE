function [S_uv,T_uv,R_uv,B_uv,h_uv,eps,a] = VinfThe2B_plane(vInf,B_theta,mu)
%--------------------------------------------------------------------------
%   Computes the geometry of the B-plane using the vertex angle 'B_theta',
%   the hyperbolic excess velocity 'vInf' and the gravitational parameter
%   'mu'.
%
%   Returns the B-plane coordinate system unit vectors 'S_uv', 'T_uv' and
%   'R_uv', as well as the unit vector of the B-vector 'B_uv', the angular
%   momentum unit vector 'h_uv', the energy of the trajectory 'eps' and its
%   semi-major axis 'a'.
%--------------------------------------------------------------------------
%   Form:
%   [S_uv,T_uv,R_uv,B_uv,h_uv,eps,a] = VinfThe2B_plane(vInf,B_theta,mu)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   vInf        (3,:)    km/s       hyperbolic excess velocity
%   B_theta     (1,:)    rad        vertex angle of B-vector and T-unit
%                                   vector
%   mu          (1,:)    km^3/s^2   gravitational parameter of central body
%
%   ------
%   Output
%   ------
%   S_uv        (3,:)    -          unit vector S
%   T_uv        (3,:)    -          unit vector T-axis of B-plane
%   R_uv        (3,:)    -          unit vector R-axis of B-plane
%   h_uv        (3,:)    -          angular momentum unit vector
%   eps         (1,:)    km^2/2^2   energy of trajectory
%   a           (1,:)    km         semi-major axis 
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 13.01.2020 |  A. Probst    | First revision
%*************************************************************************%

% S_uv, along direction of hyperbolic entry velocity
S_uv = vInf./vecnorm(vInf,2,1);

% T_uv, in xy-plane of body fixed reference frame
[T_uv] = T_unitVector(S_uv);

% R_uv, to complete coordinate system of B-Plane RTS, origin: center of
% body
[R_uv] = R_unitVector(S_uv,T_uv);

% rotation matrix from B-plane reference frame RTS --> planet inertial
% reference frame XYZ
rot_RTS = [R_uv T_uv S_uv];

% B-vector, in RST / b-plane reference frame
B_RTS = [sin(B_theta);cos(B_theta);0];

% B-vector unit vector
B_uv = rot_RTS * B_RTS;

% angular momentum unit vector
h_uv = cross(B_uv,S_uv,1);

% energy of trajectory eps, km^2/s^2
eps = 1/2 * vecnorm(vInf,2,1)^2;

% semi-major axis a, km
a = (-1) * mu/(2 * eps);

% % semi-major axis a, km
% a = (-1) * mu/vecnorm(vInf,2,1)^2;

