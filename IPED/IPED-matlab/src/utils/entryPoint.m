function [altitude,npoint,elts,states,idx,diff] = entryPoint(elts,idx,...
                                            radii,hEntry,n,et,diff)
%--------------------------------------------------------------------------
%   Defines the trajectory point at entry with an entry altitude 'hEntry' 
%   using 'idx' (altitude > hEntry) and 'idx+1' (altitude < hEntry).
%
%   Returns index 'idx' for entry point with an accuracy 'diff' < 0.1 km,
%   and the respective altitudes 'altitude' and surface points 'npoints'
%   of the trajectory arc represented by conical elements 'elts' and 
%   cartesian states 'states'.
%--------------------------------------------------------------------------
%   Form:
%   [altitude,npoint,elts,states,idx,diff] = entryPoint(elts,idx,radii, ...
%                                            hEntry,n,et,diff)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   elts        (8,:)    km,rad     SPICE conical elements of traj.
%                                   state(s)
%   idx         (3,1)    -          index of closest point of trajectory 
%                                   to entry altitude
%   radii       (:,3)    km         planet's radius at equator
%   hEntry      (1,1)    km          entry altitude above planet's surface
%   n           (1,1)    -          number of states to be propagated until
%                                   periapsis
%   et          (1,1)    sec        Ephemeris time in sec after J2000
%   diff        (1,1)    km         optional: difference of hEntry and  
%                                   trajectory altitude at idx
%
%   ------
%   Output
%   ------
%   altitude    (1,1)    km         altitude above planet's surface of 
%                                   point idx
%   npoint      (3,1)    km         planet's surface point of point idx
%   elts        (8,:)    km,rad     SPICE conical elements of traj.
%                                   state(s)
%   states      (6,:)    km,km/s    State vector of propagated trajectory
%   idx         (3,1)    -          index of state vector of trajectory 
%                                   point at entry / closest appraoch
%   diff        (1,1)    km         difference of hEntry and trajectory 
%                                   altitude at idx
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 18.03.2020 |  A. Probst    | First revision
% v1.1    | 10.11.2021 |  A. Probst    | Adding and updating comments
%*************************************************************************%

if nargin == 6
    diff = 10;
end

while diff > 0.1

    %% Iterate between state idx and idx+1
    % linear distribution 
    flag = false;
    % computing trajectory states
    [states,elts] = ...
        computingOrbStates(elts(:,idx),n,flag,et,elts(6,idx),elts(6,idx+1));
    
    %% Defining point of trajectory with closest altitude to entry altitude
    
    % Defining trajectory radius points
    r = states(1:3,:);
    
    % Calculating nearest points on ellipsoid to states and altitudes above
    % ellipsoid
    [npoint,altitude] = cspice_nearpt(r,radii(1),radii(2),radii(3));
   
    % defining index of last point before entry
    idx = find(altitude - hEntry < 0,1) - 1;
    
    % difference in altitude and entry altitude
    diff = altitude(1,idx) - hEntry;
    
end

end

