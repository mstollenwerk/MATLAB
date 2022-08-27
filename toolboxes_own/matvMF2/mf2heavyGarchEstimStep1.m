function [eparamStep1, nLogLStep1] = mf2heavyGarchEstimStep1(r,RC,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[p,~,T] = size(RC);

options = optimoptions('fminunc', 'Display', 'final', 'MaxFunEval', 1e4);

%% Step 1
m_max = 63;
nLogLStep1 = NaN(m_max,1);
eparamStep1 = NaN(m_max,p+4);

if nargin ==3
    x0Step1 = varargin{1}(1:end-2);
else
    x0Step1 = [0.05*ones(1,p), 0.2, 0.6, 0.2, 0.6]';
end
    
parfor m = 1:m_max
    obj_funStep1 = @(param) mf2heavyGarchLikeRecStep1(param,r,RC,m);

    eparamStep1(m,:) = fminunc(obj_funStep1,x0Step1,options);
    
    nLogLStep1(m) = obj_funStep1(eparamStep1(m,:)');
end