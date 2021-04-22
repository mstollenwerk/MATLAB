function [ nLogL, logLcontr, Sigma_, S, varargout ] = ...
    gas_ncF_likeRec_staticOmega( param, p, q, data_RC, mhg_precision )

[k,~,T] = size(data_RC);
k_ = k*(k+1)/2;
%% Parameter Transformation   
idx = [
    1;
    1 + k_;
    1 + k_ + k^2*p;
    1 + k_ + k^2*p + k^2*q;
    1 + k_ + k^2*p + k^2*q + k_;
    1 + k_ + k^2*p + k^2*q + k_ + 1;	
    1 + k_ + k^2*p + k^2*q + k_ + 1 + 1;
];       

chol_gas_intrcpt = param( idx(1) : idx(2) - 1 );
gas_arch_matrix  = param( idx(2) : idx(3) - 1 );
gas_garch_matrix = param( idx(3) : idx(4) - 1 );
Omega_           = param( idx(4) : idx(5) - 1 );
df_F_1           = param( idx(5) : idx(6) - 1 ); 
df_F_2           = param( idx(6) : idx(7) - 1 );

gas_intrcpt = ivechchol(chol_gas_intrcpt);
gas_arch_matrix = reshape(gas_arch_matrix,k,k,p);
gas_garch_matrix = reshape(gas_garch_matrix,k,k,q);
Omega_ = ivechchol(Omega_);

if nargout == 5

	param_out.chol_gas_intrcpt = chol_gas_intrcpt;
	param_out.gas_intrcpt = gas_intrcpt;
	param_out.gas_arch_matrix = gas_arch_matrix;
	param_out.gas_garch_matrix = gas_garch_matrix;
	param_out.Omega_ = Omega_;
	param_out.df_F_1 = df_F_1;
	param_out.df_F_2 = df_F_2;    
	param_out.all = param;

	varargout{1} = param_out;

end
	
%% Data Storage
Sigma_ = NaN(k,k,T+1);
S = NaN(k,k,T);
%% Recursion
ini_Sigma = mean(data_RC,3);
logLcontr = NaN(T,1);
for tt=1:T+1
    Sigma_(:,:,tt) = gas_intrcpt;
    if ~all(gas_arch_matrix(:)==0)
        for jj = 1:p
            if (tt-jj) <= 0
                % Left in the code for completeness:
                Sigma_(:,:,tt) = Sigma_(:,:,tt) + gas_arch_matrix(:,:,jj)*zeros(k)*gas_arch_matrix(:,:,jj)'; 
            else
                Sigma_(:,:,tt) = Sigma_(:,:,tt) + gas_arch_matrix(:,:,jj)*S(:,:,tt-jj)*gas_arch_matrix(:,:,jj)';
            end
        end
    end
    if ~all(gas_garch_matrix(:)==0)
        for jj = 1:q
            if (tt-jj) <= 0
                Sigma_(:,:,tt) = Sigma_(:,:,tt) + gas_garch_matrix(:,:,jj)*ini_Sigma*gas_garch_matrix(:,:,jj)';
            else
                Sigma_(:,:,tt) = Sigma_(:,:,tt) + gas_garch_matrix(:,:,jj)*Sigma_(:,:,tt-jj)*gas_garch_matrix(:,:,jj)';
            end
        end   
    end
    if tt < T+1
        
        % Log-Likelihoods and Scores
        try
            [~, logLcontr(tt), ScoreF] = matvncFBOPSlike( ...
                Sigma_(:,:,tt), ...
                Omega_, ...
                df_F_1, ...
                df_F_2, ...            
                mhg_precision, ...
                data_RC(:,:,tt) ...
            );
        catch
            nLogL = inf;
            logLcontr = NaN;
            return
        end            
        % Scaling the Scores as in Opschoor
		% Opschoor et al. take derivative ignogring that Sigma_ is symmetric, whereas matvFscore takes symmetry into account.
        Z = .5*(ones(k) + eye(k)) .* ivech(ScoreF.Sigma_); 
        
        S(:,:,tt) =  2/(df_F_1 + 1) * Sigma_(:,:,tt) * Z * Sigma_(:,:,tt);

    end
end
%% Likelihood
nLogL = -sum(logLcontr);
end