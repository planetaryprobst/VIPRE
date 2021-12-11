function testCase_Uranus(spicepath)

%--------------------------------------------------------------------------
%  Test case: Uranus
%
%  Computes entry conditions for an arrival at the planet Uranus using an
%  interplanetary trajectory.
%   
%  Requires the path to the generic SPICE kernel folder as input.
%
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   spicepath   str                 path to generic SPICE kernel folderr
%
%   ------
%   Output
%   ------
%   plots
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 06.12.2021 |  A. Probst    | First release
%*************************************************************************%

if( nargin == 0 )
  disp('Define path to the generic SPICE kernel folder.')
end

if ischar(spicepath)
    IPED(spicepath,'Uranus')
else
    disp('The input path needs to be a string.')
end
