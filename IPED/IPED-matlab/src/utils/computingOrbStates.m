function [states,elts] = computingOrbStates(elts0,n,flag,et,M0,MX)
%--------------------------------------------------------------------------
%   Computes orbital states based on one set of conical elements elts0 
%   and the number of following steps n. 
%
%   Returns the orbital states 'states' including the respective conical 
%   elements elts between mean anomaly of -30 deg and periapsis 0 deg, if
%   not otherwise specified. 
%--------------------------------------------------------------------------
%   Form:
%   [states,elts] = computingOrbStates(elts0,n,flag,et,M0,MX)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   elts0       (8,:)    km,rad     SPICE conical elements of traj. state
%   n           (1,1)    -          number of states to be propagated until
%                                   periapsis
%   flag        boolean  -          True for logarithmic decreasing
%                                   interval steps when approaching 
%                                   periapsis, false for equal sized steps
%   et          (1,1)    sec        Ephemeris time in sec after J2000
%   M0          (1,1)    rad        optional, default for [] is -30 rad
%   MX          (1,1)    rad        optional, default for [] is 0 rad
%
%   ------
%   Output
%   ------
%   states      (6,:)    km,km/s    State vector of propagated trajectory
%   elts        (8,:)    km,rad     SPICE conical elements of traj. states
%                                   as propagated
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 12.03.2020 |  A. Probst    | First revision
% v1.1    | 10.11.2021 |  A. Probst    | Renaming and commenting
%*************************************************************************%

if nargin < 5
    % mean anomaly of first state, in rad
    M0 = -30;
    % mean anomaly of last state, in rad
    MX = 0;
else
    % do nothing
end

% Duplicating set of conic elements of state @ periapsis 
elts = repmat(elts0,1,n);

if flag
    % Varying mean anomaly in unequal steps, with decreasing step size
    % approaching MX starting at M0
    M = (M0 - MX) * (logspace(log10(1.5),log10(0.5),n) - 0.5) + MX;
else
    % Varying mean anomaly in equal steps between M0 and MX
    M = (M0 - MX) * [1:-1/(n-1):0] + MX;
end

% Updating states with new mean anomaly with MX as the last state
elts(6,:) = M;

% RV states of orbit until MX
states = cspice_conics(elts,repmat(et,1,size(elts,2)));

end





