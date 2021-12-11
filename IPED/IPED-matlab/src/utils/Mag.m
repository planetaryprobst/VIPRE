function m = Mag( u )
%--------------------------------------------------------------------------
%   Computes the magnitudes 'm' of each column of a given 3-by-n matrix 'u'
%   with each column representing a vector.
%
%   Returns a row vector of the magnitudes of each column. 
%--------------------------------------------------------------------------
%   Form:
%   m = Mag( u )
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   u            (:,:)  Vectors
%
%   -------
%   Outputs
%   -------
%   m            (:)   Corresponding magnitudes
%
%*************************************************************************%
% Language: MATLAB R2019b (OSX)
% Author: Alena Probst
% History:
% Version |    Date    |     Name      | Change history
% v1.0    | 10.11.2021 |  A. Probst    | First revision
%*************************************************************************%

m = sqrt(sum(u.^2));

