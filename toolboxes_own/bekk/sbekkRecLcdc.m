function [ S, fcst ] = sbekkRecLcdc( param, P, Q, data )
%SBEKKREC returns positive definite matrices from BEKK recursion

n = size(data{1},1);
T = max(size(data));
%% Transformations
L = reshape(param(1:n^2),n,n);
d = param(n^2+1:n^2+n);
a = param(n^2+n+1:n^2+n+P);
b = param(n^2+n+P+1:n^2+n+P+Q);
Sc = L*diag(d)*L';
intrcpt=(1-sum(a)-sum(b))*Sc;
%% Initialization
S = cell(1,T+1);
ini = mean(cat(3,data{:}),3);
%% Recursion
for t=1:T+1
    S_ = intrcpt;
    for p=1:P
        if t-p<1
            S_ = S_ + a(p)*ini;
        else
            S_ = S_ + a(p)*data{t-p};
        end
    end
    for q=1:Q
        if t-q<1
            S_ = S_ + b(q)*ini;
        else
            S_ = S_ + b(q)*S{t-q};
        end
    end
    S{t} = S_;
end
fcst = S{end};
S = S(1:end-1);
end
