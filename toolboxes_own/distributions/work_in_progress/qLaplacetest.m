clear
clc

df = 3;
Sigma_ = [0.5 1.5; 2 0.5];
Sigma_ = Sigma_ * Sigma_';

for ii = 1:10
    [ qLaplace_rnd ] = qLaplacernd(ones(1e4,1)*df, repmat(Sigma_,1,1,1e4));
    [ eparam{ii}, tstats{ii}, logL{ii}, optimoutput{ii} ] = qLaplaceest( qLaplace_rnd, [], 'Display', 'iter-detailed', 'PlotFcn', 'optimplotx' );
end

for ii = 1:10
    eSigma(:,:,ii) = eparam{ii}.Sigma_;
    e_df(ii) = eparam{ii}.df;
end

figure
boxplot(squeeze(eSigma(1,1,:)))
title('(1,1)')

figure
boxplot(squeeze(eSigma(2,1,:)))
title('(2,1)')

figure
boxplot(squeeze(eSigma(2,2,:)))
title('(2,2)')

figure
boxplot(e_df)
title('df')