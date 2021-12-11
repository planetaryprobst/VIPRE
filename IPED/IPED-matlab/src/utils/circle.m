function [x,y,z] = circle(o, r, n)
%--------------------------------------------------------------------------
%   Rotates a circle with radius 'r', center coordinates 'o' from a plane
%   with normal vector 'n0' = [0,0,1] into a plane with normal vector 'n'.   
%   
%   Returns the coordinates 'x','y','z' of the points of the rotated 
%   circle.
%--------------------------------------------------------------------------
%   Form:
%   [x,y,z] = circle(o, r, n)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   o           (3,1)     km        position vector of circle origin
%   r           (1,1)     km        radius 
%   n           (3,1)     km        vector normal to plane
%
%   ------
%   Output
%   ------
%   x           (:,1)      km       x-coordinates of new circle
%   y           (:,1)      km       y-coordinates of new circle
%   z           (:,1)      km       z-coordinates of new circle
%
%*************************************************************************%
% Language: MATLAB R2017b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 18.04.2018 |  A. Probst    | First revision
% v1.1    | 10.11.2021 |  A. Probst    | Updating comments
%*************************************************************************%

% points of a circle in original plane
t = linspace(0,2*pi);
x = r * cos(t);
y = r * sin(t);
z = 0*t;
pnts = [x;y;z];

% unit normal vector of original plane
n0 = [0;0;1]; 

% unit normal vector of final plane
n = n/vecnorm(n);

% theta is the angle between the two normal vectors
c = dot(n0,n) / ( vecnorm(n0)*vecnorm(n) ); % cos(theta)
s = sqrt(1-c*c);                        % sin(theta)
u = cross(n0,n) / ( vecnorm(n0)*vecnorm(n) ); % rotation axis as unit vector
%u = u/vecnorm(u) % ... as unit vector

C = 1-c;

% rotation matrix
ROT = [u(1)^2*C+c, u(1)*u(2)*C-u(3)*s, u(1)*u(3)*C+u(2)*s
    u(2)*u(1)*C+u(3)*s, u(2)^2*C+c, u(2)*u(3)*C-u(1)*s
    u(3)*u(1)*C-u(2)*s, u(3)*u(2)*C+u(1)*s, u(3)^2*C+c];


% Rotated points
newPnts = o + ROT*pnts;

%% OUTPUT
x = newPnts(1,:)';
y = newPnts(2,:)';
z = newPnts(3,:)';

