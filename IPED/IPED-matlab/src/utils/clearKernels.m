function success = clearKernels()
%--------------------------------------------------------------------------
%   Clears previously loaded SPICE kernels for IPED-matlab computations. 
%
%   Returns a success boolean.
%--------------------------------------------------------------------------
%   Form:
%   success = clearKernels()
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%     -   
%
%   ------
%   Output
%   ------
%   success     bool     -          Boolean, true for successful execution
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 20.05.2020 |  A. Probst    | First revision
%*************************************************************************%

%% clear spice kernels

cspice_kclear

success = true;


end

