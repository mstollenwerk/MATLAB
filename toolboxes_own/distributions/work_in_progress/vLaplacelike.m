function [ nLogL, logLcontr ] = vLaplacelike( data_mat, Sigma_ )
%VLAPLACELIKE
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
%      [1] Kotz, Kozubowski and Podgorski (2001), p. 250, 6.5.3
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 20.02.2020
%
% DEPENDENCIES:
%
%%
[p,~,n] = size(Sigma_);
[n_,p_] = size(data_mat);
if n ~= n_ || p ~= p_
    error('Dimensions of data_mat and Sigma_ do not agree.')
end

nu_ = (2-p)/2;
logLcontr = NaN(n,1);
for ii = 1:n
    c = log(2) - p/2*log(2*pi) - .5*log(det(Sigma_(:,:,ii)));
    Q = data_mat(ii,:)*(Sigma_(:,:,ii)\data_mat(ii,:)');
    log_kernel = nu_/2*log(Q/2) + log(besselk(nu_, sqrt(2*Q)));
          
    logLcontr(ii) = c + log_kernel;
end
nLogL = -sum(logLcontr);
end
