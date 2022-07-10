function R = matvsitWishrnd( Sigma_, n, nu, N )
%MATVTWISHRND
%
% USAGE:
%   
%
% INPUTS:
%   
%
% OUTPUTS:
%   
%  See also 
%
% COMMENTS:
%   
% REFERENCES:
%      [1]                    
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 21.03.2021
%
% DEPENDENCIES:

Omega_ = matvStandardize('itWish',Sigma_,[n; nu]);
R = matvitWishrnd( Omega_, n, nu, N );

end