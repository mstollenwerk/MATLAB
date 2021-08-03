function [ nLogL, score ] = mvstBOPSlike( Sigma_, skew, df, x)
%MVSTBOPSLIKE Negative log-likelihood for the multivariate skew-t as 
%   defined in Bodnar, Okhrin, Parolya and Stollenwerk (2020) [BOPS].
%
% USAGE:
%   [ nLogL, score ] = mvstBOPSlike( Sigma_, alpha_, nu_, x)
%
% INPUTS:
%   SIGMA_  - Array (p by p). Symmetric p.d. parameter matrix, 
%             regulates covariance. 
%   SKEW    - Array (p by 1). Parameter vector, regulates skewness.
%   DF      - Double. Degrees of freedom parameter.
%   X       - Array (p by 1). Data.
%
% OUTPUTS:
%   NLOGL   - Double. Negative log-likelihood value.
%   SCORE   - Struct. Fields as parameter names. Contains derivatives
%             of log-likelihood (matvFlike) w.r.t. parameters.
%
% COMMENTS:
%   This distribution is slightly different from Azzalini (2013).
%
% REFERENCES:
%   [1] Bodnar, Okhrin, Parolya and Stollenwerk (2020) [BOPS]S
%   [2] Azzalini, Adelchi. The skew-normal and related families. 
%           Vol. 3. Cambridge University Press, 2013.
%
% DEPENDENCIES:
%
%  See also
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 04.06.2020

%% Input Checking
narginchk(4,5);
if size(skew,2) > size(skew,1)
    skew = skew';
end
if size(x,2) > size(x,1)
    x = x';
end
p = size(x,1);
%% Log-likelihood computation
invSig = inv(Sigma_);
Q = x'*invSig*x;
kern_ = 1 + Q / (df-2);

term1 = gammaln((df+p)/2) - gammaln(df/2) - p/2*log((df-2)*pi);
term2 = -log(det(Sigma_))/2;
term3 = -(df+p)/2*log(kern_);
mvtpart = term1 + term2 + term3;

nLogL = - mvtpart;

% Boundary case of central t distribution.
if any(skew ~= 0)
    t_arg = skew'*x*sqrt(df + p)/sqrt(df - 2 + Q);
    % This code for the tcdf is taken from the matlab function tcdf [R2020,
    % Statistics and Machine Learning Toolbox]
    % They quote:
    %   References:
    %      [1] M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
    %      Functions", Government Printing Office, 1964, 26.7.
    %      [2] L. Devroye, "Non-Uniform Random Variate Generation",
    %      Springer-Verlag, 1986
    %      [3] E. Kreyszig, "Introductory Mathematical Statistics",
    %      John Wiley, 1970, Section 10.3, pages 144-146.
    xsq = t_arg^2;
    tcdf_y = betainc(xsq ./ (df + p + xsq), 0.5, (df + p)/2, 'upper') / 2;
    % For x > 0, F(x) = 1 - F(-|x|).
    if t_arg > 0
        tcdf_y = 1 - tcdf_y;
    end

    tcdfpart = log( tcdf_y );

    nLogL = -log(2) - mvtpart - tcdfpart;
end
%% Score computation
if nargout == 2
    
    % Boundary case of central t distribution.
    if any(skew ~= 0)
        c_nu = sqrt(df + p) / sqrt(df - 2);

        tpdf_y = ...
            exp(gammaln((df + p + 1) / 2) - gammaln((df + p)/2)) ./ ...
            (sqrt((df + p)*pi) .* ...
            (1 + (t_arg .^ 2) ./ (df + p)) .^ ((df + p + 1)/2));

        w = c_nu / kern_ * ...
            (c_nu + tpdf_y / tcdf_y * skew' * x / (df-2) /sqrt(kern_) );

        dev_ = x*x'*w - Sigma_;

        scoreSig = .5*invSig*dev_*invSig;
        
        score.skew = tpdf_y/tcdf_y*sqrt(df + p)/sqrt(df - 2 + Q)*x;
    else % stupidly adjusted:
        w = (df + p)/(df - 2)/kern_;

        scoreSig = .5.*(invSig*(x*x')*invSig.*w - invSig);
        score.skew = zeros(p,1);
    end
    
    score.Sigma_ = 2*scoreSig - scoreSig.*eye(p); %double off diagonals to take symmetry of Sigma into account in the derivative.

end

end

