function [x,y,z] = plot_planetSphere(radii,nrFacets)
%--------------------------------------------------------------------------
%   Creates and ellipsoid in the planet's shape. Takes the number of facets
%   'nFacets' and the planet's dimensions in x,y,z in 'radii' as input.
%
%   Returns the ellipsoid point matrices in 'x','y' and 'z'.
%
%--------------------------------------------------------------------------
%   Form:
%   [x,y,z] = plot_planetSphere(radii,nrFacets)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   radii      (3,:)    km        planet's radius in x,y,z
%   nFacets    (1)      -         number of facets of ellipsoid
%
%   ------
%   Output
%   ------
%   vRot       (:,1)    km/s       rotational speed at latitude
%   R_lat      (:,1)    km         radius at latitude lat 
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 12.01.2020 |  A. Probst    | First revision
%*************************************************************************%
 
% ellipsoidal references frame axes, body fixed
u = [1,0,0;0,1,0;0,0,1];

% ellipsoid point matrices
[x,y,z]=Ellipsoid(radii,u,nrFacets);

end