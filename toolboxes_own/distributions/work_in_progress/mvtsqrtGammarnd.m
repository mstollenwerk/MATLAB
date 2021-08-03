function [ X ] = mvtsqrtGammarnd( Sigma_, df_t, lambda_, N )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 27.07.2020
p = size(Sigma_,1);
%% Random Draws
mvn = mvnrnd(zeros(p,1), eye(p), N);
gam1 = gamrnd(2*df_t, 2/df_t, N, 1);
gam2 = gamrnd(lambda_, 1/lambda_, N, 1);

X = mvn./repmat(sqrt(gam1),1,p).*repmat(sqrt(gam2),1,p);
