function wxy = Cross( w, y )
%--------------------------------------------------------------------------
%   Computes the vector cross product. The number of columns of 'w' and 'y'
%   can be both > 1 and equal.
%
%   One can have one column and the other any number of columns.
%--------------------------------------------------------------------------
%   Form:
%   wxy = Cross( w, y )
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   w                 (3)    Vector
%   y                 (3)    Vector
%
%   -------
%   Outputs
%   -------
%   wxy               (3)    Vector cross product of w and y
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 25.05.2020 |  A. Probst    | First revision
%*************************************************************************%


cW = size(w,2);
cY = size(y,2);

if( cW == cY )
  wxy = [w(2,:).*y(3,:)-w(3,:).*y(2,:);...
         w(3,:).*y(1,:)-w(1,:).*y(3,:);...
	   w(1,:).*y(2,:)-w(2,:).*y(1,:)]; 
		 
elseif( cW == 1)
  wxy = [w(2)*y(3,:)-w(3)*y(2,:);...
         w(3)*y(1,:)-w(1)*y(3,:);...
	   w(1)*y(2,:)-w(2)*y(1,:)]; 
		 
elseif( cY == 1)
  wxy = [w(2,:)*y(3)-w(3,:)*y(2);...
         w(3,:)*y(1)-w(1,:)*y(3);...
	   w(1,:)*y(2)-w(2,:)*y(1)]; 
else
  error('w and y cannot have different numbers of columns unless one has only one column')
end

