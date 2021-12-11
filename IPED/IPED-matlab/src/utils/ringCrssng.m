function [flag] = ringCrssng(rings,state)
%--------------------------------------------------------------------------
%   Checks if input state 'state' lies within the radial distance of the
%   rings of the planet, with the geometric information of the planet given
%   in 'rings'. Assumes a circular ring formation (particles of rings on
%   circular orbits).
%
%   Returns a 'flag', 'true' for the state located within the rings,
%   'false' for the state located outside of the rings.
%--------------------------------------------------------------------------
%   Form:
%   [flag] = ringCrssng(rings,state)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   rings       container km        mapping container that includes the
%                                   inner and outer borders of the planet's
%                                   rings
%   states      (6,:)    km,km/s    State vector of input trajectory
%
%   ------
%   Output
%   ------
%   flag        bool     -          true: trajectory state lies within one
%                                   of the rings
%                                   false: trajectory state lies outside of
%                                   the rings
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 29.04.2020 |  A. Probst    | First revision
%*************************************************************************%

% first assumption: state in equatorial plane is NOT located in one of the 
% rings
flag = false;

% magnitude of radius vector of state = distance to center of planet
% state(3,1) == 0
r = sqrt(state(1,1)^2 + state(2,1)^2);

% extracting keys from map container
keySet = keys(rings);

% for all rings in the set
for i = 1:length(keySet)
    
    borders = rings(keySet{i});
    
    if and(borders(1)<r,r<borders(2))
        flag = true;
    end
    
end

end     