clear
clc

df_1 = 3;
df_2 = 5;
Sigma_ = [0.5 1.5; 2 0.5];
Sigma_ = Sigma_ * Sigma_';

for ii = 1:10
    [ qLaplace_rnd(:,:,:,ii) ] = ...
        matvtWishrnd(Sigma_, df_1, df_2, 1e4);
    [ eparam{ii}, tstats{ii}, logL{ii}, optimoutput{ii} ] = ...
        matvtWishest( qLaplace_rnd(:,:,:,ii), [], 'Display', 'iter-detailed' );
end

for ii = 1:10
    eSigma(:,:,ii) = eparam{ii}.Sigma_;
    e_df_1(ii) = eparam{ii}.df_1;
    e_df_2(ii) = eparam{ii}.df_2;
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
boxplot(e_df_1)
title('df_1')

figure
boxplot(e_df_2)
title('df_2')