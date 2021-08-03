function A = matviWishrnd( Sigma_, df, N )
%MATVIWISHRND
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
% DEPENDENCIES:
%
% Michael Stollenwerk
% michael.stollenwerk@live.com

p = size(Sigma_,1);
A = NaN(p,p,N);
for ii=1:N

    A(:,:,ii) = iwishrnd( Sigma_, df );
    
end

end

