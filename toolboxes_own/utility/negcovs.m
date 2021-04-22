function [sumnegcovs, relnegcovs, relnegcovstotal] = negcovs(data_mat)
%NEGCOVS counts the number of negative covariances in covariance matrix
%time sereis
%   Detailed explanation goes here

[k,~,n] = size(data_mat);

sumnegcovs = NaN(n,1);
relnegcovs = NaN(n,1);
for ii = 1:n
    dta = data_mat(:,:,ii);
    sumnegcovs(ii) = sum(dta(tril(true(k),-1))<0);
    relnegcovs(ii) = mean(dta(tril(true(k),-1))<0);
end
relnegcovstotal = sum(sumnegcovs)/(k*(k-1)/2)/n;
end 