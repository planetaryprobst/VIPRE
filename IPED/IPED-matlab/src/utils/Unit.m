function u = Unit( v )
%--------------------------------------------------------------------------
%   Unitizes vectors by column.
%--------------------------------------------------------------------------
%   Form:
%   u = Unit( v )
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   v            (:,n)     Vectors
%
%   -------
%   Outputs
%   -------
%   u            (:,n)     Unit vectors
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 12.03.2020 |  A. Probst    | First revision
%*************************************************************************%

[n,p]  = size(v);

if( n == 1 )
  m = v;
else
  m = sqrt(sum(v.^2));
end

% start with NaN
u = NaN(n,length(m));

% check that magnitudes are not zero
k = find( m > 0 );

if( ~isempty(k) )
  if (p==1)
    % Single column
    u = v/m;
  else
    for j = 1:n
      u(j,k) = v(j,k)./m(k);	
    end 
  end
end

