function phi = RV2FPA( r, v, mu )
%--------------------------------------------------------------------------
%   Computes the flight path angle 'phi' for a cartesian state given with
%   its position vector 'r', its velocity vector 'v' and the gravitational
%   parameter 'mu'.
%--------------------------------------------------------------------------
%   Form:
%   phi = FlightPathAngle( r, v, mu )
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   r               (6,:) [r;v] or (3,1) r
%   v               (3,:) v
%   mu              (1,1) Gravitational parameter (only needed to plot)
%
%   -------
%   Outputs
%   -------
%   phi             (1,:) Flight path angle
%
%--------------------------------------------------------------------------
%   Reference: Vallado, D. A. Fundamentals of Astrodynamics and Applications
%              McGraw-Hill. pp. 138-139.
%--------------------------------------------------------------------------


if( nargin < 1 )
  [r,v] = RVOrbGen( [24472 0 0 0 0.5 0] );
  FlightPathAngle(r,v);
  return;
end

if( nargin == 1 )
  v = r(4:6,:);
  r = r(1:3,:);
end

if( nargin < 3 )
    mu = 3.98600436e5;
end

rMag   = Mag(r);
vMag   = Mag(v);

% Must deal with the zero velocity case
%--------------------------------------
j      = find(vMag == 0);
v(:,j) = Unit(r(:,j))*1e-6;
vMag(j) = Mag(v(:,j));

%phi = pi/2 - acos( Dot(r,v) ./ (Mag(r).*Mag(v)) );

rDot   = sum( r.*v )./rMag;

h      = Mag(cross(r,v));

cosPhi = h./(rMag.*vMag);
sinPhi = rDot./vMag;
phi    = atan2( sinPhi, cosPhi );

if( nargout == 0 )
  disp('here')
  n  = length(phi);
  nu = zeros(1,n);
  for k = 1:n
      el    = RV2El(r(:,k), v(:,k), mu );
      nu(k) = M2Nu(el(5),el(6));
  end
  k = find(nu<0);
  nu(k) = nu(k) + 2*pi;
  Plot2D(nu*180/pi,phi*180/pi,'Sample','Phi (deg)','Flight Path Angle')
  clear phi;
end

