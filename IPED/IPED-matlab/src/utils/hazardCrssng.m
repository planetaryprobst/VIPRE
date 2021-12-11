function [state,safe] = hazardCrssng(states,elts,et,radii,rings)
%--------------------------------------------------------------------------
%   Checks if the trajectory is passing through a hazard.
%    
%   Returns a boolean 'safe' and a state vector 'state' that represents the
%   closest trajectory state in that leg to the equatorial plane. 
%
%   Based on a trajectory state set 'states' and its equivalent set of
%   conical elements 'elts', it is checked if the trajectory leg in 
%   'states'/'elts' crosses the equator at a location where there are no
%   rings.
% 
%--------------------------------------------------------------------------
%   Form:
%   [state,safe] = hazardCrssng(states,elts,et,radii,rings)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   states      (6,:)    km,km/s    State vector of input trajectory
%   elts        (8,:)    km,rad     SPICE conical elements of traj. states
%                                   as computed
%   et          (1,1)    sec        Ephemeris time in sec after J2000
%   radii       (3,1)    km         radii of planet in X, Y, Z direction
%   rings       container km        mapping container that includes the
%                                   inner and outer borders of the planet's
%                                   rings
%
%   ------
%   Output
%   ------
%   state       (6,:)    km,km/s    State vector of computed trajectory
%                                   of position that is closest to 
%                                   equatorial plane
%   safe        bool     -          true:  trajectory is safe
%                                   false: trajectory is unsafe
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 29.04.2020 |  A. Probst    | First revision
%*************************************************************************%


% boolean indicating trajectory is safe 
safe = true;

%% Equatorial Plane Crossing

% Check if trajectory crosses the equatorial plane in the given trajectory
% set, flag = true for yes, flag = false for no 
[state,~,flag] = eqtCrssng(states,elts,et);

% if trajectory does cross the equatorial plane in this set
if flag
    % check if the crossing happens inside the planet
    if and(flag,(state(1,1)/radii(1))^2 + (state(2,1)/radii(2))^2 <=1) 
        % then trajectory is considered as safe but the state shall be
        % ignored in the output
        state(1:6,:) = NaN(6,1); 
    else
        % check if the trajectory state lies within one of the rings
        if ringCrssng(rings,state)
            % then trajectory is considered as unsafe
            safe = false;
        end
    end
else
    % trajectory is safe
    % do nothing
end

end
 