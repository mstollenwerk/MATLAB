function [ S, fcst ] = sbekkRec( params, covMdata, P, Q )
%SBEKKREC returns positive definite matrices from BEKK recursion
%   [ S ] = sbekkRec( params, data, p, q ) requires input parameters,
%   PARAMS, to be the non-zero entries of the lower cholesky factor of the
%   intercept matrix followed by the coefficients for the lagged DATA,
%   followed by the coefficients for the lagged recursions:
%   [vech(chol(Sc,'lower'))', a_1,...,a_p, b_1,...,b_q]. DATA must be
%   input as cell of length T, where the t'th cell entry contains the t'th
%   (realized) covariance matrix (for quicker recursion evaluation). P are
%   the number of lagged DATA, Q are the number of lagged RECURSION.
%
%   [ S, fcst ] = sbekkRec( params, data, p, q ) also returns the one
%   period ahead forecast.

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 29.08.2016

[N,~,T] = size(covMdata);
%% Transformations
C = ivech(params(1:N*(N+1)/2),'lower');
C = C*C';
if any(eig(C,'vector'))<=0
    error('The intercept matrix is not positive definite');
end
a = params(N*(N+1)/2 + 1 : N*(N+1)/2+P);
b = params(N*(N+1)/2+P+1 : end);
%% Initialization
initial = mean(covMdata,3);
S = cell(1,T+1);
%% Recursion
for t=1:T+1
    S_ = C;
    for p=1:P
        if t-p<1
            S_ = S_ + a(p)*initial;
        else
            S_ = S_ + a(p)*covMdata(:,:,t-p);
        end
    end
    for q=1:Q
        if t-q<1
            S_ = S_ + b(q)*initial;
        else
            S_ = S_ + b(q)*S{t-q};
        end
    end
    S{t} = S_;
end
fcst = S{end};
S = S(1:end-1);
end