function [eparam,optimoutput] = fmincon_Rieszperm(p,obj_fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,varargin)
%FMINCON_RIESZPERM implements Algorithm 11 of Blasques, Lucas, Opschoor and
%Rossini (2021) - Tail Heterogeneity for Dynamic Covariance-Matrix-Valued
%Random Variables - The F-Riesz Distribution.

rng(1); % Random Seed
n_perm = 4;
perm_ = NaN(n_perm,p);
perm_(1,:) = 1:p;
for ii = 2:2*n_perm % Create matrix of originally supplied and 2*n_perm other permutations. 2*n_perm, because there is a chance that randperm generates the same permuation again in the loop.
    perm_(ii,:) = randperm(p);
end
perm_ = unique(perm_, 'rows', 'stable'); % Remove duplicate perm_s
perm_ = perm_(1:n_perm,:); % Now only take the first n_perm rows.
opt_like = NaN(n_perm,1);

parfor ini_perm = 1:n_perm
    
    perm_try = perm_(ini_perm,:); % This has to be done for parfor loop to work.
    [eparam_,optimoutput_] = ...
            my_fmincon(...
                @(x) obj_fun(x, perm_try), ...
                x0, ...
                A,b,Aeq,beq,lb,ub,nonlcon,varargin{:} ...
            );

    for ii = 1:p

        i_star = find(perm_(ini_perm,:)==ii);
        perm1_ii = perm1(perm_(ini_perm,:),ii);    

        lls = NaN(p,1);
        for jj = 1:p
            lls(jj) = obj_fun(eparam_,perm1_ii(jj,:));
        end
        i_opt = find(lls == min(lls));

        if i_opt ~= i_star
            perm_(ini_perm,:) = perm1_ii(i_opt,:);
            
            perm_try = perm_(ini_perm,:); % This has to be done for parfor loop to work.
            [eparam_,optimoutput_] = ...
                    my_fmincon(...
                        @(x) obj_fun(x, perm_try), ...
                        eparam_, ...
                        A,b,Aeq,beq,lb,ub,nonlcon,varargin{:} ...
                    );
        end

    end
    
    optimoutput_.perm_ = perm_(ini_perm,:);
    eparam{ini_perm} = eparam_;
    optimoutput{ini_perm} = optimoutput_;
    opt_like(ini_perm) = optimoutput_.history.fval(end);
    
end

opt_id = find(opt_like==min(opt_like));
eparam = eparam{opt_id(1)};
optimoutput = optimoutput{opt_id(1)};

end