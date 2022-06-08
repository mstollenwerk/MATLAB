function R = matviWishrnd( Omega_, df, N )
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

p = size(Omega_,1);
C = chol(Omega_,'lower');

R = NaN(p,p,N);
for ii=1:N

    BU = BarlettU(ones(p,1)*df);
    L = C/BU';
    R(:,:,ii) = L*L';
    
end

end

