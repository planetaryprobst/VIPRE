function [omega] = planetRotation(planet)
%--------------------------------------------------------------------------
%   Computes the rotational period of the planet 'omega' in rad/sec.
%--------------------------------------------------------------------------
%   Form:
%   [omega] = planetRotation(planet)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   planet      str      -          SPICE code or string for planet
%
%   ------
%   Output
%   ------
%   omega       (1,1)    rad/s      rotational period of the planet
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 12.03.2020 |  A. Probst    | First revision
%*************************************************************************%

% prime meridian in deg (element 1) and deg/day (element 2)
w = cspice_bodvrd(planet,'PM',3); 

% rotational period of the planet in rad/sec
omega = w(2) * cspice_rpd() / cspice_spd();


