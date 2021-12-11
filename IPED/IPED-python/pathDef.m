function success = pathDef(homePath)
%--------------------------------------------------------------------------
%   Defines the necessary directories and paths needed for the IPED
%   computations. Adds the paths to the Matlab path.
%
%   Uses the home directory path to the IPED directory as input.
%
%   Returns a success boolean.
%--------------------------------------------------------------------------
%   Form:
%   IPED(spicepath,body)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   homePath    str    path     home directory path to IPED directory:
%                               /home/user/IPED/
%
%   ------
%   Output
%   ------
%   
%   success     boolean         true for successful completion
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 07.12.2021 |  A. Probst    | First release
%*************************************************************************%
    
    % defining paths for Matlab 
    addpath(genpath([homePath,'NAIF/']))
    addpath(genpath([homePath,'IPED-matlab/']))
    
    success = true;

end