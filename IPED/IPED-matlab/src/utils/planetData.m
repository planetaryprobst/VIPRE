function [mu,radii,f,rings] = planetData(body)
%--------------------------------------------------------------------------
%   Load the planetary data of the body specified in the input string
%   'body'.
%
%   Returns the gravitational parameter 'mu', the body's 'radii' for the x,
%   y, and z axis of the body, the flattening factor 'f' and the ring
%   geometry in 'rings'.
%--------------------------------------------------------------------------
%   Form:
%   [mu,radii,f,rings] = planetData(body)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   planet      str      -          SPICE code or string for body
%
%   ------
%   Output
%   ------
%   mu          (1,1)    km^2/s^2   Gravitaional parameter
%   radii       (3,1)    km         raddi for x, y and z axis of body
%   f           (3,1)    -          vectorized flattening factor
%   rings       map      km         ring geometry in km
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 25.05.2020 |  A. Probst    | First revision
%*************************************************************************%


%% Body Data

% gravitational parameter
mu = cspice_bodvrd( body, 'GM', 1 );

% ellipsoidal radii (1) = x, (2) = y, (3) = z
radii = cspice_bodvrd( body, 'RADII', 3 );

% flattening factor in vectorization / -
f = [1;1;radii(3)/radii(1)]; 


%% Rings

if strcmp(body,'Saturn') || strcmp(body,'Uranus') || strcmp(body,'Neptune')

    % ring definition
    rings = true;

end

if rings 
    
    % Creating key-calue set for planet's rings, in km
        
    if strcmp(body,'Saturn')
        % reference radius for ring dimensions, km
        r = max(radii(1),radii(2));
        % radii and dimensions of Saturnian rings in km
        keySet = {'D','C','B','A','F','G','E'};
        valueSet = {[1.11 1.236] * r [1.239 1.526] * r [1.526 1.95] * r ...
            [2.03 2.27] * r [2.32 2.32] * r [2.754 2.874] * r...
            [2.987 7.964] * r};
        rings = containers.Map(keySet,valueSet);
        
    else
        % for Uranus and Neptune, the total euqatorial plane is assumed as
        % a 
        keySet = {'0'};
        valueSet = {[min(radii(1),radii(2)),8 * max(radii(1),radii(2))]};
        rings = containers.Map(keySet,valueSet);
          
    end
end
    


