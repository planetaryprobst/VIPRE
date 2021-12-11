function [vRot,R_lat] = rotationalSpeed(lat,omega,radii)
%--------------------------------------------------------------------------
%   Computes the rotational speed of a planet, given its 'radii' in 
%   x,y,z directions and a rotational velocity 'omega' at a given
%   latitude 'lat'.
%
%   Returns the rotational speed 'vRot' of a point on the surface at
%   latitude 'R_lat'. Returns a vector if 'lat' is a row vector.
%   
%--------------------------------------------------------------------------
%   Form:
%   [vRot,R_lat] = rotationalSpeed(lat,omega,radii)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   lat        (:,1)    rad       latitude of point on planet's surface
%   omega      (1,1)    rad/s     rotational speed of planet
%   radii      (3,:)    km        planet's radius at equator
%
%   ------
%   Output
%   ------
%   v_rot      (:,1)    km/s       rotational speed at latitude
%   R_lat      (:,1)    km         radius at latitude lat 
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 12.01.2020 |  A. Probst    | First revision
% v1.1    | 10.02.2020 |  A. Probst    | Change input unit for radius from
% m to km
%*************************************************************************%
 
% radius at latitude lat and flattening factor f
[R_lat,~] = radiusNonSphericalPlanet(radii,lat);

% rotational speed of body at latitude lat
vRot  = omega .* R_lat .* cos(lat);

end