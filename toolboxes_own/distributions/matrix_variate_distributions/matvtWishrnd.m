function R = matvtWishrnd( Omega_, n, nu, T )
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
%   This is the distribution of the quadratic form XX', where the columns
%   of the p by df_n matrix X follow multivariate central t distributions
%   with the same degree of freedom df_t and the same dispersion matrix 
%   Sigma_.
%   
% REFERENCES:
%      [1]                    
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.08.2020
%
% DEPENDENCIES:

p = size(Omega_,1);
R = NaN(p,p,T);
for ii=1:T

    W = wishrnd( Omega_, n );
    chi2_ = chi2rnd(nu); 
    R(:,:,ii) = W/(chi2_/nu);
    
    % Alternatively:
%     gamma_ = gamrnd( df_t/2, 2/df_t );
%     A(:,:,ii) = W./gamma_;
    
end

end

