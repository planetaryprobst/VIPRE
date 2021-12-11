function success = loadKernels(spicepath)
%--------------------------------------------------------------------------
%   Loads necessary SPICE kernels for IPED-matlab computations. Takes the
%   path to the generic SPICE folder as input.
%
%   Returns a success boolean.
%--------------------------------------------------------------------------
%   Form:
%   success = loadKernels(spicepath)
%--------------------------------------------------------------------------
%
%   -----
%   Input
%   -----
%   spicepath   str                 Path to SPICE generic_kernels directory
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


%% LOAD STANDARD SPICE KERNELS

if nargin == 0
    spicepath = ['../NAIF/naif.jpl.nasa.gov/pub/naif/generic_kernels/'];
end

% Planetary Ephemeris
file = 'spk/planets/de432s.bsp';
cspice_furnsh( [spicepath,file] )

% Planetary GM Data
file = 'pck/gm_de431.tpc';
cspice_furnsh( [spicepath,file] ) 

% Planet Orientation and Radii
file = 'pck/pck00010.tpc';
cspice_furnsh( [spicepath,file] )

% Leapseconds
file = 'lsk/naif0012.tls';
cspice_furnsh( [spicepath,file] )

success = true;

end

