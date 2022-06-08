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

p = size(Sigma_,1);
Omega_ = (n-p-1)*Sigma_;
R = matvitWishrnd( Omega_, n, nu, N );

end