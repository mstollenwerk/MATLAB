function A = matvitWishrnd( Omega_, n, nu, N )
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

p = size(Omega_,1);
A = NaN(p,p,N);
for ii=1:N

    W = iwishrnd( Omega_, nu );
    chi2_ = chi2rnd(n); 
    A(:,:,ii) = W*(chi2_/n);
    
    % Alternatively:
%     gamma_ = gamrnd( df_t/2, 2/df_t );
%     A(:,:,ii) = W.*gamma_;
    
end

end

