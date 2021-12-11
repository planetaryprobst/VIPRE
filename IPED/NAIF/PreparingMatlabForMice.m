addpath('/Documents/NAIF/mice/src/mice/')
addpath('/Documents/NAIF/mice/lib/' )

%% FUNCTIONS OF INTEREST

% % Converting r,v in long/lat 
% 
% filepath = 'naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc';
% 
%       %
%       % load PCK file
%       %
%       cspice_furnsh( filepath )
% 
%       %
%       % retrieve Mars radii
%       %
%       radii = cspice_bodvrd( 'MARS', 'RADII', 3 );
%       f = (radii(1) - radii(3))/radii(1);
% 
%       %
%       % package 3D vector
%       %
%       %vec = [3373.850; -351.034; -117.267]
%         
%       vec =  [-299761;-440886;-308712]
%       %
%       % compute planetocentric coordinates
%       %
%       [ pcr, pclon, pclat ] = cspice_reclat( vec )
% 
%       %
%       % compute planetodetic coordinates
%       %
%       [ pdlon, pdlat, pdalt ] = cspice_recgeo( vec, radii(1), f ) 
% 
%       %
%       % compute planetographic coordinates
%       %
%       [ pglon, pglat, pgalt ] = cspice_recpgr( 'MARS', vec, radii(1), f )

Example(1):

      %
      % Determine the osculating elements of the moon wrt the
      % Earth at some arbitrary time in the J2000 inertial frame.
      %
      % Load the meta kernel listing the needed SPK, PCK, LSK
      % kernels.
      %
      cspice_furnsh('standard.tm' )

      %
      % Convert the time string to ephemeris time
      %
      et = cspice_str2et( 'Dec 25, 2007' );

      %
      % Make the cspice_spkezr call to retrieve the state of the
      % moon wrt the Earth in J2000.
      %
      [state, ltime] = cspice_spkezr( 'Moon', et, 'J2000', 'LT+S', 'EARTH' );

      %
      % cspice_oscelt requires body mass information, so load a
      % mass PCK kernel.
      %
      cspice_furnsh( '/kernels/gen/pck/masses3.tpc' )

      %
      % Read the gravitational parameter for Earth.
      %
      mu = cspice_bodvrd( 'EARTH', 'GM', 1 );

      %
      % make the cspice_oscelt call to convert the state 6-vector
      % to the elts 8-vector. Note: the  cspice_bodvrd returns
      % data as arrays, so to access the gravitational parameter
      % (the only value in the array), we use mu(1).
      %
      elts = cspice_oscelt( state, et, mu(1) );

      %
      % Output the elts vector in a column format.
      %
      txt = sprintf( '%24.8f\n', elts );
      disp( txt)

      %
      % It's always good form to unload kernels after use,
      % particularly in MATLAB due to data persistence.
      %