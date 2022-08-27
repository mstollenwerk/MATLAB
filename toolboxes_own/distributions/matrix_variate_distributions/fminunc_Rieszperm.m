function [eparam,optimoutput] = fminunc_Rieszperm(p,obj_fun,x0,varargin)
%FMINCON_RIESZPERM implements Algorithm 11 of Blasques, Lucas, Opschoor and
%Rossini (2021) - Tail Heterogeneity for Dynamic Covariance-Matrix-Valued
%Random Variables - The F-Riesz Distribution.

tic

rng(1); % Random Seed
if p <= 10
    if p == 2
        n_perm = 2;
    elseif p == 3
        n_perm = 3;
    else
        n_perm = 4;
    end
else % It gets too computationally expensive to have multiple starting points.
    n_perm = 1;
end
perm_ = NaN(n_perm,p);
% If initial permutation is given, use this one.
for ii = 1:numel(varargin)
    if isstring(varargin{ii}) || ischar(varargin{ii})
        if strcmp(varargin{ii},'perm')
            perm_(1,:) = varargin{ii+1};
            perm_input_idx = ii;
        end
    end
end
if exist('perm_input_idx','var')
    varargin = varargin([1:perm_input_idx-1,perm_input_idx+2:numel(varargin)]);
end
% If initial permutation is not given use 1:p.
if isnan(perm_(1,1))
    perm_(1,:) = 1:p;
end
for ii = 2:2*n_perm % Create matrix of originally supplied and 2*n_perm other permutations. 2*n_perm, because there is a chance that randperm generates the same permuation again in the loop.
    perm_(ii,:) = randperm(p);
    while isinf(obj_fun(x0, perm_(ii,:))) || isnan(obj_fun(x0, perm_(ii,:)))
        warning('Random Permutation didnt work with x0.')
        perm_(ii,:) = randperm(p);
    end
end
perm_ = unique(perm_, 'rows', 'stable'); % Remove duplicate perm_s
perm_ = perm_(1:n_perm,:); % Now only take the first n_perm rows.

opt_like = NaN(n_perm,1);

parfor ini_perm = 1:n_perm
    
    perm_try = perm_(ini_perm,:); % This has to be done for parfor loop to work.
    [eparam_,optimoutput_] = ...
            my_fminunc(...
                @(x) obj_fun(x, perm_try), ...
                x0, ...
                varargin{:} ...
            );

    for ii = 1:p

        i_star = find(perm_(ini_perm,:)==ii);
        perm1_ii = perm1(perm_(ini_perm,:),ii);    

        lls = NaN(p,1);
        for jj = 1:p
            lls(jj) = obj_fun(eparam_,perm1_ii(jj,:));
        end
        i_opt = find(lls == min(lls));
        i_opt = i_opt(1);
        
        if i_opt ~= i_star
            perm_(ini_perm,:) = perm1_ii(i_opt,:);
            
            perm_try = perm_(ini_perm,:); % This has to be done for parfor loop to work.
            try
                [eparam_,optimoutput_] = ...
                        my_fminunc(...
                            @(x) obj_fun(x, perm_try), ...
                            eparam_, ...
                            varargin{:} ...
                        );
            catch
                [eparam_,optimoutput_] = ...
                        my_fminunc(...
                            @(x) obj_fun(x, perm_try), ...
                            x0, ...
                            varargin{:} ...
                        );
            end

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

optimoutput.estimation_time_overall = toc;

end