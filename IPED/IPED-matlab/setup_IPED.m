%--------------------------------------------------------------------------
%   IPED - Impact of the location of the Planetary Entry probe on
%   spacecraft and mission Design
%--------------------------------------------------------------------------
%   Script to compute entry conditions for all safe entry
%   opportunities of a planetary entry probe being released from an
%   interplanetary trajectory. 
%   
%   The entry conditions are computed for all safe entries available for a
%   chosen interplanetary trajectory. A trajectory is considered safe if it
%   doesn't cross a ring area.
%   
%   It creates several plots to visualize the results.
%
%   IPED is part of the software package VIPRE, consisting of VAPRE and 
%   IPED,  a software to Visualize the Impact of the PRobe Entry location 
%   on the science, mission adn spacecraft design.
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 10.11.2021 |  A. Probst    | First release
%*************************************************************************%

clearvars
close all


%% USER INPUT

% individual path to NAIF folder
mypath = '/Users/aprobst/Documents/';

% path to the generic spice kernel folder
spicepath = [mypath,'NAIF/naif.jpl.nasa.gov/pub/naif/',...
    'generic_kernels/'];

% body name: choose a body available in SPICE
body = 'Saturn';

%% FUNCTION CALL
IPED(spicepath,body)
