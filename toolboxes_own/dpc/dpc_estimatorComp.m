function [ eparamComp, e_d, fcstd ] = dpc_estimatorComp( covMdata, diagLRL, compSpec, P, Q, x0Comp )
%DPCCAWESTIMATORCOMP estimates components recursion in DPCCAW
%   DPCCAWESTIMATORCOMP 

[N,~,T] = size(covMdata);
%% Defaults
if isempty(Q)
     Q=1;
end
if isempty(P)
    P=1;
end
%% Storage and Reshape x0 for individual deployment
e_d = NaN(N,T);
fcstd=NaN(N,1);
x0Comp = reshape(x0Comp,N,[]); % Each x0Comp row now corresponds to starting values for one component.
edc = sort(eig(mean(covMdata,3)),'descend');
%% Minimizations
% options = optimoptions('fminunc','Display','off','Algorithm','quasi-newton');
options = optimoptions('fmincon','Display','off');
if strcmpi(compSpec,'garch') || isempty(compSpec)
    %x0
    if isempty(x0Comp)
        x0Comp=ones(100,P+Q)*.01;
    end        
    %Restrictions
    AAgarch = [-eye(P+Q)  % Non-negativity constraint on eigenvalues to ensure p.d. RC matrix.
            ones(1,Q+P)]; % Stationarity 
    bbgarch = [zeros(1,P+Q) 1-eps];
    %Minimization
    eparamComp = NaN(size(x0Comp));
    for i=1:N
        eparamComp(i,:) = fmincon(@(params)cawRVqlike(garchRec([(1-sum(params))*edc(i) params],diagLRL(i,:),P,Q),diagLRL(i,:)), x0Comp(i,:), AAgarch, bbgarch,[],[],[],[],[],options);
%         eparamComp(i,:) = fminunc(@(params)cawRVqlike(garchRec([(1-sum(params))*edc(i) params],diagLRL(i,:),P,Q),diagLRL(i,:)), x0Comp(i,:), options);
        [e_d(i,:), fcstd(i)] = garchRec([(1-sum(eparamComp(i,:)))*edc(i) eparamComp(i,:)],diagLRL(i,:),P,Q);
    end
elseif strcmpi(compSpec,'har')
    %x0
    if isempty(x0Comp)
        x0Comp=ones(100,4)*.1;
    end  
    %Restrictions
    AAhar = [-eye(4)  % Non-negativity constraint on eigenvalues to ensure p.d. RC matrix.
            ones(1,4)]; % Stationarity 
    bbhar = [zeros(1,4) 1-eps];
    %Minimization
    eparamComp = NaN(size(x0Comp));
    for i=1:N
        eparamComp(i,:) = fmincon(@(params)cawRVqlike(har1dRec([(1-sum(params))*edc(i) params],diagLRL(i,:)),diagLRL(i,:)), x0Comp(i,:),AAhar,bbhar,[],[],[],[],[],options);
        [e_d(i,:), fcstd(i)] = har1dRec([(1-sum(eparamComp(i,:)))*edc(i) eparamComp(i,:)],diagLRL(i,:));
    end
% elseif strcmpi(compSpec,'midas20m12lag')
%     if ~isempty(x0Comp) % b/c midasRV allows for empty x0 input.
%         x0Comp = reshape(x0Comp,N,5); % Each x0Comp row now corresponds to starting values for one component.
%     end
%     edc = sort(eig(mean(covMdata,3)),'descend');
%     e_d = [NaN(100,20*12) repmat(edc,1,T-20*12)]; % Storage and Initialization. Set data used for midas filter initialization with NaN values. Initialize to unconditional mean edc. This ensures all components that are not included in the Factor-DPCCAW are set to the proper value.
%     eparamComp = NaN(N,5); % Storage
% %     tic
%     for i=1:compNfac
%         if ~isempty(x0Comp) %UGLY, BUT TOO LAZY NOW.
%             [eparamComp(i,:), ~, e_d(i,:), ~, ~ ] = midasRV(diagLRL(i,:),12,20,x0Comp(i,:));
%         else
%             [eparamComp(i,:), ~, e_d(i,:), ~, ~ ] = midasRV(diagLRL(i,:),12,20,[]);
%         end
%     end
% %     disp(strcat('Components parameters maximum likelihood estimation finished after ', num2str(toc), ' seconds.'));
%     eparamComp = eparamComp(:)';
else error('Only GARCH and HAR specification supported so far.')
end
eparamComp = eparamComp(:)';
end