function score = matvncFBOPSscore_rnd(df_F_1, df_F_2, Sigma_, Omega_, RC)
    p = size(Sigma_,1);
    
    mcrep = 1e5;
    N = 100;
    
    c = df_F_1/(df_F_2 - p - 1);
    inv_Sigma_new = Sigma_/RC*Sigma_ + c*Sigma_;
    Sigma_new = inv(inv_Sigma_new);
    Omega_new = c*(inv_Sigma_new\Omega_)/inv_Sigma_new;
    invRC = inv(RC);
    invSigma = inv(Sigma_);
    Theta_new = Sigma_new\Omega_new;
    
    mean_nom = NaN(p,p,N);
    parfor ii = 1:N
        
        matvncWish_rnd = matvncWishrnd( ...
            df_F_1, ...
            Sigma_new, ...
            Theta_new, ...
            mcrep ...
        );

        nom = NaN(p,p,mcrep);
        for jj = 1:mcrep

            V = matvncWish_rnd(:,:,jj);
                    
            Z_V = .5*( V*Sigma_*invRC + invRC*Sigma_*V + c*V );
            detV = det(V);

            nom(:,:,jj) = Z_V * detV^(df_F_2/2);

        end
        mean_nom(:,:,ii) = mean(nom,3);
        
    end
    
    mean_detV_power_h = ...
        matvncWish_Edet(df_F_1, Sigma_new, Theta_new, df_F_2/2);
    Z = (.5* df_F_1 + df_F_2)*invSigma ...
        + .5* invSigma*Omega_*invSigma ...
        - mean(mean_nom,3)./mean_detV_power_h;
    
    score.Sigma_ = 2*Z - Z.*eye(p);
    
end