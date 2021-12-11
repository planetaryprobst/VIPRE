function [R,f] = radiusNonSphericalPlanet(radii,lat)
%--------------------------------------------------------------------------
%   Returns the planet radius 'R' and ellipticity 'f' at a geocentric
%   latitude 'lat' for a non-spherical planet (oblate in z) with an equal
%   equatorial radius in x and y direction, given with the vector 'radii'.
%--------------------------------------------------------------------------
%   Form:
%   [R,f] = radiusNonSphericalPlanet(radii,lat)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   radii       (3,1)    km     SPICE planet radii [rX;rY;rZ]
%   lat         (:,1)    rad    geocentric latitude in radians
%
%   ------
%   Output
%   ------
%   R           (:,1)     km     radius at latitude delta
%   f           (1,1)     -      ellipticity of the planet / flattening
%   factor
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 12.01.2020 |  A. Probst    | First revision
%*************************************************************************%
   

% ellipticity of the planet
f = 1 - radii(3)/radii(1); 

% Radius at geocentric latitude delta, km
% derived from equation of an ellipsoid: (a) 1 = x^2/radii(1)^2 +
% z^2/radii(3)^2, (b) x = R * cos(lat) und (c) z = R * sin(lat)
R = sqrt(radii(1)^2 * radii(3)^2 ./ ...
        ((radii(3) .* cos(lat)).^2 + (radii(1) .* sin(lat)).^2));

end
