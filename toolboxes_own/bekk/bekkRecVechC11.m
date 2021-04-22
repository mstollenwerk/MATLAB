function [ S, fcst ] = bekkRecVechC11( param, data )
%SBEKKREC returns positive definite matrices from BEKK recursion

[n,~] = size(data{1});
T=max(size(data));
%% Transformations
C = ivech(param(1:n*(n+1)/2),'sym'); % Passing in vech from NOT chol!
if any(eig(C,'vector'))<=0
    error('The intercept matrix is not positive definite');
end
A = reshape(param(n*(n+1)/2 + 1 : n*(n+1)/2+n^2),n,n); % Actually you could only input the vech form since they are sym pd.
B = reshape(param(n*(n+1)/2+n^2+1 : end),n,n);
%% Initialization
S = cell(1,T+1);
ini = mean(cat(3,data{:}),3);  
S{1}=C + B*ini*B + A*ini*A;
%% Recursion
for t=2:T+1
    S{t} = C + B*S{t-1}*B + A*data{t-1}*A;
end
fcst = S{end};
S = S(1:end-1);
end
