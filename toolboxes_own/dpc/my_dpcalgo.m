function [ eparam, logQL, eS, eL, e_d ] = my_dpcalgo(covMdata, cP, cQ, bl_l)
[n,~,T]=size(covMdata);
covMdata_cell = mat2cell(covMdata(:,:),n,ones(1,T)*n);
%% dpcestimator to find starting point
tic
[ eparam, fval_, ~, ~, ~, ~, ~, ~, ~ ] = dpc_estimator( covMdata, [], 1, 1, [], cP, cQ, [] );
disp('===================================================================')
disp('')
disp('DPC_ESTIMATOR')
toc
disp('')
disp('===================================================================')
fval=fval_-1; % to start the while loop
run = 0;
while fval_-fval>0.01
    run=run+1;
    fval = fval_;
    [eLc, edc, rest_eparam] = dpcpartrans(n,eparam);
    %% Estiamate dc with restrictions in groups of bl_l
    if bl_l == n
        AA = [-eye(bl_l-1) zeros(bl_l-1,1)]+[zeros(bl_l-1,1) eye(bl_l-1)];
        bb = zeros(1,bl_l-1);
        lb = zeros(bl_l,1);
        % Optimization
        format long
        edc = fmincon(@(edc)sdpcLikeRecLcdc( [eLc(:)' edc rest_eparam], covMdata_cell,1,1,cP,cQ ), edc',AA,bb,[],[],lb,[],[],optimoptions('fmincon','Display','iter-detailed','FunctionTolerance',1e-2,'StepTolerance',1e-3,'MaxIterations',50,'UseParallel','always'));
        disp('===================================================================')
        disp('')
        disp(strcat('Eigenvalues all'))
        toc
        disp('')
        disp('===================================================================')
        format short
    else        
        for i=1:bl_l:n
            tic
            if i==1 % Restrictions
                AA_ = [-eye(bl_l-1) zeros(bl_l-1,1)]+[zeros(bl_l-1,1) eye(bl_l-1)];
                AA = [zeros(1,bl_l-1) -1; AA_];
                bb = [-edc(i+bl_l) zeros(1,bl_l-1)];
            elseif i==(n-bl_l-1)
                AA_ = [-eye(bl_l-1) zeros(bl_l-1,1)]+[zeros(bl_l-1,1) eye(bl_l-1)];
                AA = [1 zeros(1,bl_l-1); AA_];
                bb = [edc(i-1) zeros(1,bl_l-1)];
            else 
                AA_ = [-eye(bl_l-1) zeros(bl_l-1,1)]+[zeros(bl_l-1,1) eye(bl_l-1)];
                AA = [zeros(1,bl_l-1) -1; 1 zeros(1,bl_l-1); AA_];
                bb = [-edc(i+bl_l) edc(i-1) zeros(1,bl_l-1)];
            end
            % Optimizations
            format long
            edc(i:i+bl_l-1)=fmincon(@(edc_i)sdpcLikeRecLcdc( [eLc(:)' edc(1:i-1)' edc_i edc(i+bl_l:end)' rest_eparam], covMdata_cell,1,1,cP,cQ ), edc(i:i+bl_l-1)',AA,bb,[],[],[],[],[],optimoptions('fmincon','Display','iter-detailed','FunctionTolerance',1e-2,'StepTolerance',1e-3,'MaxIterations',50,'UseParallel','always'));
            disp('===================================================================')
            disp('')
            disp(strcat('Eigenvalues dc ',num2str(i),'to',num2str(i+bl_l-1)))
            toc
            disp('')
            disp('===================================================================')
            format short
        end
    end
    %% Estimate lA and lB
    tic
    rest_eparam(1:2)=fminunc(@(x)sdpcLikeRecLcdc( [eLc(:)' edc x rest_eparam(3:end)], covMdata_cell,1,1,cP,cQ ), rest_eparam(1:2),optimset('Display','iter-detailed','UseParallel','always'));
    disp('===================================================================')
    disp('')
    disp('A and B of Loading Recursion')
    toc
    disp('')
    disp('===================================================================')
    %% Estimate cA and cB
    tic
    eparam_ = dpcpartrans(eLc, edc, rest_eparam);
    eparamComp = dpc_estComp( covMdata, eparam_(1:n*(n+1)/2+2), 1, 1, cP, cQ, eparam_(n*(n+1)/2+3:end));
    eparam_ = [eparam_(1:n*(n+1)/2+2) eparamComp];
    disp('===================================================================')
    disp('')
    disp('As and Bs of Component Recursion')
    toc
    disp('')
    disp('===================================================================')
    %% New fval
    fval_ = sdpcLikeRec( eparam_, covMdata, cP, cQ );
    fval_ = -fval_;
    disp('===================================================================')
    disp('')
    disp('Fval diff:')
    disp(fval_)
    disp(fval)
    disp(fval_-fval)
    disp('End run')
    disp(run)
    disp('')
    disp('===================================================================')
    if fval_-fval<=0.05
       break
    elseif 0.05 < fval_-fval
       eparam = eparam_;
    end
end
%% Outputs
[ nLogL, ~, eS, eL, e_d, ~ ] = sdpcLikeRec( eparam, covMdata, cP, cQ );
logQL = -nLogL;
end