function [state,dist,flag] = eqtCrssng(states,elts,et)
%--------------------------------------------------------------------------
%   Searches the trajectory state with the minimum distance to the 
%   equatorial plane for a given trajectory state set 'states' and its 
%   equivalent set of conical elements 'elts' by iterating and computing
%   states between the two states that are enclosing the crossing. 
%   
%   Returns the trajectory state 'state' that lies in the equatorial plane
%   with a maximum deviation of 0.1 km in z-direction. If no crossing
%   appears in the given trajectory set, the state with the closest 
%   distance is returned. A boolean flag 'flag' is 'true' for the first, 
%   'false' for the latter case.
% 
%--------------------------------------------------------------------------
%   Form:
%   [state,dist,flag] = eqtCrssng(states,elts,et)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   states      (6,:)    km,km/s    State vector of input trajectory
%   elts        (8,:)    km,rad     SPICE conical elements of traj. states
%                                   as computed
%   et          (1,1)    sec        Ephemeris time in sec after J2000
%
%   ------
%   Output
%   ------
%   state       (6,:)    km,km/s    State vector of new computed trajectory
%                                   states that is closest to equatorial 
%                                   plane
%   dist        (1,1)    km         vertical distance of point on
%                                   equatorial plane (NaN for false cases) 
%   flag        bool     -          true:  state in equatorial plane 
%                                   false: state not in equatorial plane
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 20.04.2020 |  A. Probst    | First revision
%*************************************************************************%

% size of given input trajectory set
setSize = size(states,2);

% number of iterations between state p and state k
n = 10;

% flag is true for trajectory sets that cross the equatorial plane
flag = true;

%% Determine Equatorial Plane Crossing

% find closest distance to equator plane
[dist,p] = min(abs(0-states(3,:)));

% set flag if closest point of trajectory is first point in states
if p == 1
    first = true;
else
    first = false;
end

% determining crossing point of trajectory with an accuracy of 0.1 km 

% while the distance to the equatorial plane is above 0.1 km
while dist > 0.1
    % iterate between trajectory state below and above the equatorial
    % plane until the threshold is undercut
    
    % if the closest point is the first point in state set
    if p == 1
        % increase it by one
        p = 2;
    end
    
    % if the closest trajectory state and the trajectory state before are 
    % both on the same hemisphere of the planet (+Z or -Z)
    if sign(states(3,p)) == sign(states(3,p-1))
        % and the closest trajectory state is the last point in states 
        if p == size(states,2)
            % and if it is within the first iteration
            if p == setSize
                % then the last point is periapsis 
                % and the equator crossing happens after periapsis and
                % therefore after entry
                k = NaN;
            else
                % else there has been a numerical error. To fix that, the
                % last value of the iteration set is ignored
                k = 1;
                p = p - 1;
            end
        % else if it is the first point in states within the first
        % iteration
        elseif first
            % then all other point thereafter are further away from
            % the equator than the first one and there is no
            % equator crossing
            p = 1;
            k = NaN;
        else    
            % the equator crossing happens before periapsis, between p and 
            % p+1
            k = p + 1;
        end
    else
        % then the equator crossing happens between p and p-1
        k = p - 1;
    end
    % if equator crossing happens after entry or closest point is first
    % time in states within the first iteration
    if isnan(k)
        % no ring crossing is possible before entry
        dist  = NaN;
        flag = false;

    else
        % find the state when trajectory crosses the equatorial plane by
        % computing points between p and k
        [states,elts] = ...
            computingOrbStates(elts(:,min(p,k)),n,false,et, ...
            elts(6,min(p,k)),elts(6,max(p,k)));
        % and find the state with closest distance to equator plane
        [dist,p] = min(abs(0-states(3,:)));
    end
end

% if distance criterion meets z-axis component
if abs(states(3,p)) <= 0.1
    % define the equator crossing state when dist <= 0.1 km
    state = states(:,p);
else
    % for all other states set position vector to NaN
%     state = [NaN;NaN;NaN;states(4:6,p)];
    state = NaN(6,1);
end

end


            