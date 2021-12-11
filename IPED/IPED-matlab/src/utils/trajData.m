function [trajectories,vInf,epoch] = trajData(body)
%--------------------------------------------------------------------------
%   Loads the trajectory data of the planet specified in the input string
%   planet from a Matlab data file with format '.mat'.
%
%   Returns a table 'trajectories' with the trajectory information, the
%   hyperbolic entry velocity for all trajectories in the vector 'vInf',
%   and a vector 'epoch' with the arrival time for each trajectory.
%
%--------------------------------------------------------------------------
%   Form:
%   [trajectories,vInf,epoch] = trajData(planet)
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
%   trajectories table   SI units   interplanetary trajectories for body
%   vInf        (3,:)    km/s       hyperbolic entry veclocity per
%                                   trajectory
%   epoch       (1,:)    days       arrival time
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 25.05.2020 |  A. Probst    | First revision
%*************************************************************************%

% NAIF ID associated to body
[ ~, found] = cspice_bodn2c( body );

if found==0
    disp('The input planet is not a valid planet!')
    return
end

%% Trajectory Data

try 
    % Interplanetary Trajectories
    trajectories = load(fullfile('data', [body,'.mat'])).(body);
catch ME
    % error message 
    disp('The trajectory data is missing!')
    return
end

% hyperbolic incoming velocity
vInf = [trajectories{:,3} trajectories{:,4} trajectories{:,5}]';

epoch = trajectories{:,2}';
        
end

    


