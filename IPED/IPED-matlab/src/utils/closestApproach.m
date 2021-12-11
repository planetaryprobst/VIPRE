function [altitude,npoint,idx,flag] = closestApproach(r,radii,hEntry)
%--------------------------------------------------------------------------
%   Computes the point of closest approach of a trajectory with respect to
%   the surface of a planet.
%
%   Returns (1) the 'altitude' above the surface of a planet with 'radii' 
%   X,Y,Z for all trajectory points in 'r', (2) the coordinates 'npoint' of
%   the projections of each trajectory point in 'r' on the planet's 
%   surface, (3) the index 'idx' of 'r' for the closest approach of the
%   trajectory to the planet, and (4) a boolean 'flag' indicating if r(idx)
%   is below the defined entry altitude 'hEntry', with 'true' indicating 
%   an entry trajectory, 'false' indicating a flyby trajectory.
%--------------------------------------------------------------------------
%   Form:
%   [altitude,npoint,idx,flag] = closestApproach(r,radii,hEntry)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   r          (3,:)    km          matrix with trajectory position vectors
%   radii      (3,1)    km          radii of planet in X, Y, Z direction
%   hEntry     (1,1)    km          entry altitude above planet's surface
%
%   ------
%   Output
%   ------
%   altitude   (1,:)    km          altitude above planet's surface for all
%                                   trajectory points
%   npoint     (3,;)    km          planet's surface points of all
%                                   trajectory points
%   idx        (3,1)    -           index of state vector of trajectory 
%                                   point at entry / closest appraoch
%   flag       boolean  -           true for entry, false for closest
%                                   approach
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 13.03.2020 |  A. Probst    | First revision
% v1.1    | 18.03.2020 |  A. Probst    | Implementing MICE to calculate
% clossest points and altitude to ellipsoid
%*************************************************************************%

%% Defining trajectory points

% Calculating nearest points on ellipsoid to states and altitudes above
% ellipsoid
[npoint,altitude] = cspice_nearpt(r,radii(1),radii(2),radii(3));


%% Closest approach of trajectory to surface

% closest approach aka periapsis
[altitude_idx,idx] = min(altitude);

% check if closest approach / periapsis is lower than entry altitude
if altitude_idx < hEntry
    % defining index of last point before entry
    idx = find(altitude - hEntry < 0,1) - 1;
    % entry flag set to true
    flag = true;
else
    % no entry as periapsis is not below the given entry altitude 
    flag = false;
end
