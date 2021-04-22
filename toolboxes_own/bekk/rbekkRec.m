function [ S, fcst ] = rbekkRec( param, data )
%SBEKKREC returns positive definite matrices from BEKK recursion

n = size(data{1},1);
T = max(size(data));
%% Transformations
Sc = ivech(param(1:n*(n+1)/2),'lower');
Sc = Sc*Sc';
[L,d]=dpceig(Sc);

a = param(n*(n+1)/2+1:n*(n+1)/2+n);
b = param(n*(n+1)/2+n+1:n*(n+1)/2+2*n);

A = L*diag(sqrt(a))*L';
B = L*diag(sqrt(b))*L';
intrcpt=L*(diag(d)-diag(sqrt(a))*diag(d)*diag(sqrt(a))-diag(sqrt(b))*diag(d)*diag(sqrt(b)))*L';
%% Initialization
S = cell(1,T+1);
ini = mean(cat(3,data{:}),3);  
S{1}=intrcpt + B*ini*B + A*ini*A;
%% Recursion
for t=2:T+1
    S{t} = intrcpt + B*S{t-1}*B + A*data{t-1}*A;
end
fcst = S{end};
S = S(1:end-1);
end
