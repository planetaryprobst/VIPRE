function [x, y, z] = Ellipsoid( s, u, n )
%--------------------------------------------------------------------------
%   Calculates an ellipsoid given the ellipsoid semimajor axes 's', the
%   unit vectors along the axes 'u', and the number of patches 'n' (similar
%   to the sphere function from Matlab).
%
%   Returns the 'x', 'y', 'z' components of the vertices of each patch.
%--------------------------------------------------------------------------
%   Form:
%   [x, y, z] = Ellipsoid( s, u, n )
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   
%   s         (3,1)      semi-major axes, [a;b;c]
%   u         (3,3)      principal axes
%   n         (1)        number of patches  
%
%   -------
%   Outputs
%   -------
%   x         (n+1,n+1)  x-component of vertices of patch
%   y         (n+1,n+1)  y-component of vertices of patch
%   z         (n+1,n+1)  z-component of vertices of patch
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 15.05.2020 |  A. Probst    | First revision
%*************************************************************************%

if nargin == 0
  % Test the rotation of a 2D ellipsoid
  Ellipsoid( [4 2 0], Eul2Mat([0;0;-pi/6]), 20 );
  hold on
  plot3(4*sin(0:0.1:2*pi),4*cos(0:0.1:2*pi),0*[0:0.1:2*pi],'r')
  axis equal
  view(0,90);
  return;
end

% Initialize the parameter angles 
phi = (-n:2:n)/n*pi;
theta = 0.5*(-n:2:n)'/n*pi; 

costheta = cos(theta);
sintheta = sin(theta);
sinphi = sin(phi);
cosphi = cos(phi);

a = s(1);
b = s(2);
c = s(3);

x = a * costheta * sinphi;
y = b * sintheta * sinphi;
z = c * ones(n+1,1) * cosphi;

% rotates the coordinates
for i = 1:n+1
    for j = 1:n+1
      r = u*[x(i,j), y(i,j), z(i,j)]';
      x(i,j) = r(1);
      y(i,j) = r(2);
      z(i,j) = r(3);
   end
end

if nargout == 0
  NewFig('Ellipsoid')
  sH = surf(x,y,z);
  set(sH, 'EdgeColor', 'none', 'FaceColor', 'b');
  axis equal
  xlabel('X - zenith (m)')
  ylabel('Y - along-track (m)')
  zlabel('Z - cross-track (m)')
  camlight left;
  material shiny;
  lighting phong;
  view(3); 
  clear x
end
