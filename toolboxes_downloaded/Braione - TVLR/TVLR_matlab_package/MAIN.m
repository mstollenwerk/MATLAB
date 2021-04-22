%==========================================================================
%==========================================================================
%     =================================================================
%
%                 A Time-Varying Long Run HEAVY model
%
%      ================================================================
%==================        Braione (2016)     =============================
%==========================================================================
% AUTHOR: Manuela Braione, CORE-UCL, Belgium
% Contact: manuela.braione@uclouvain.be
% Last Revised: July 2016
%--------------------------------------------------------------------------
% Please cite:
% Braione, M. "A Time-Varying Long Run HEAVY model." Statistics and
% Probability Letters  (2016), pp. 36-44
%--------------------------------------------------------------------------
% Brief Code description: 
%
% This code estimates and computes direct predictions from the TVLR-HEAVY 
% model as done in Braione(2016). Full sample estimation is performed by 
% setting "forecast=0" while an out-of-sample prediction exercise is allowed by
% setting "forecast=1" with forecast horizon "h=1,5 and 10" days ahead. The
% code also generates latex tables (to visualize full-sample estimation
% results) and plots of fitted assets volatilities. All settings can be
% changed (see below).
%==========================================================================



%==========================================================================
% Load your own dataset or load the one I provide
%==========================================================================

load BJ_RC5min       %load series on realized covariance matrices
load BJ_dailyret     %load daily returns
load BJ_tickers      %load tickers if any, otherwise create one string vector called 'tickers'
load BJ_dates        %load dates (if any, otherwise see below, line 30)

% % If no dates are provided, obtain them as follows:
% ini='1-Feb-2001';     %define initial day of sample
% fin='31-Dec-2009';    %define ending day of sample
% [date,date2]=createDates(ini,fin);

Vt         = BJ_RC5min;  % realized covariance, if using another dataset change 'BJ_RC5min'
T          = size(Vt,3);  % length of sample (same for realized and daily data)
n          = size(Vt,1);  % number of assets in the dataset
DR         = BJ_dailyret; % daily returns, if different change accordingly
numlags    = 260;        % number of lags backwards, usually one year
forecast   = 1;          % if 0 then only full sample estimation is performed
outsample  = 300;        % define length of out-of-sample period
horizons   = [1 5 10];   % define multi-step prediction horizons
disp_lik   = 0;          % 1 for loglik visualization, 0 otherwise
options    = optimset('Algorithm','interior-point');
options    = optimset(options,'display','off');
latex.makedocument=1;    % for building the latex table



%==========================================================================
% Define type of empirical application to perform
%==========================================================================

if forecast==1
    % Expanding window scheme:
    % IN SAMPLE    : t= 1....(T-outsample)
    % OUT OF SAMPLE: t= (T-outsample)+1 ....T
    message    = 'OUT  OF SAMPLE EXERCISE: Expanding window scheme ';
    insample   = T-outsample;
    step       = length(horizons);
    dateOOS    = date(insample+1:end);
    computeSTD = 0;
    
else
    % Full sample estimation
    outsample  = 1;
    insample   = T;
    dateOOS    = date;
    step       = 1;
    horizons   = 1;
    computeSTD = 1;
    message    = ['FULL-SAMPLE ESTIMATION ON ' num2str(insample) ' OBSERVATIONS'];
    
end


CODE_SETTING= struct('Forecast', forecast, 'insample',insample,'outsample',...
    outsample,'n',n,'T',T,'numlags', numlags,'std_error', computeSTD,...
    'horizons',horizons,'estim_options', options,'latex_options',latex);


%==========================================================================
% Compute daily series of covariance matrices from squared returns
%==========================================================================

DRt=zeros(n,n,T);
for t=1:T
    DRt(:,:,t)=DR(t,:)'*DR(t,:);
end


disp('===================================================================')
disp('')
disp(['              ',  message  ])
disp('')
disp('===================================================================')


for H=1:step
    h=horizons(H);
    forecast_TVLR_Ht=[];
    forecast_TVLR_St=[];
    
    for K=1:outsample
        jj=K-1;
        
        if outsample>1
            disp('------------------------------------------------------');
            disp([' Estimation step ', num2str(K),' at hor = ',num2str(h),'']);
            disp('------------------------------------------------------');
        end
        
        Vin    =Vt(:,:,1:insample+jj);
        DRin   =DR(1:insample+jj,:);
        DRtin  =DRt(:,:,1:insample+jj);
        tin    =size(Vin,3);
        
        
        % Estimation
        run TVLR_HEAVY_estimation
        
        
        % Storing results at each iteration
        TVLR_estimates(K,:)=TVLR_HEAVY.stime;
        TVLR_loglik(K,:)   =TVLR_HEAVY.lik;
        TVLR_longrun(K,:)  =TVLR_HEAVY.longrun;
        TVLR_int(K,:)      =TVLR_HEAVY.int;
        TVLR_GARCH(:,:,K)  =TVLR_HEAVY.GARCH;
        TVLR_DCC(K,:)      =TVLR_HEAVY.DCC;
        TVLR_exitF(K,:)    =TVLR_HEAVY.exit;
        TVLR_time(K,:)     =TVLR_HEAVY.time;
        TVLR_aic_bic(K,:)  =[TVLR_HEAVY.aic TVLR_HEAVY.bic];
        
        if computeSTD==1
            TVLR_std(K,:)  =TVLR_HEAVY.std;
            TVLR_tstat(K,:)=TVLR_HEAVY.tstat;
            
        end
        
        % Predicting out-of-sample
        TVLR_stime = TVLR_estimates(K,:);
        
        clear mat
        [mat]=cov_aggr(Vin);
        
        if forecast==1
            disp('=============================================================')
            disp(    ' Predicting at h steps ahead')
            disp('=============================================================')
            disp(' ')
        end
        
        [TVLR_Ht,TVLR_St]= predict_TVLR_HEAVY_h(TVLR_HEAVY.stime,...
            mat,Vin,DRin,numlags,h);
        
        if K==1
            forecast_TVLR_Ht = TVLR_Ht; %predicted H_t
            forecast_TVLR_St = TVLR_St; %predicted S_t
        elseif K>1
            forecast_TVLR_Ht(:,:,end+1) = TVLR_Ht(:,:,end);
            forecast_TVLR_St(:,:,end+1) = TVLR_St(:,:,end);
        end
        
        
    end
    
    % Recover horizon of prediction
    if computeSTD==1
        TVLR_step=struct('estimates', TVLR_estimates,'loglik',...
            TVLR_loglik, 'exitF',TVLR_exitF,'aic_bic',TVLR_aic_bic,...
            'int',TVLR_int,'Longrun', TVLR_longrun,'GARCH', TVLR_GARCH,'DCC',...
            TVLR_DCC,'forecast_Ht',forecast_TVLR_Ht, 'forecast_St',...
            forecast_TVLR_St,'std',TVLR_std,'tstat',TVLR_tstat);
    else
        TVLR_step=struct('estimates', TVLR_estimates,'loglik',...
            TVLR_loglik, 'exitF',TVLR_exitF,...
            'int',TVLR_int,'Longrun', TVLR_longrun,'GARCH', TVLR_GARCH,'DCC',...
            TVLR_DCC,'forecast_Ht',forecast_TVLR_Ht, 'forecast_St',...
            forecast_TVLR_St);
    end
    
    
    %STORING: change name according to the horizon
    name=matlab.lang.makeValidName(['TVLR_',num2str(h),'step']);
    Forecast.(sprintf(['TVLR_',num2str(h),'step']))=TVLR_step;
    
    %%=========================================================================
    % Writing Latex table for full sample Likelihood values
    %%=========================================================================
    
    latexTable_general   = latex_inSample(latex,n,TVLR_step);
    
    latexTable_intercept = latex_inSample_intercept(latex,n,Vt,TVLR_step);
    
    
end

% SAVING: change name depending on application
if forecast==1
    save(['Forecast_out',num2str(outsample),'_nlags',num2str(numlags),'_hor',num2str(h),'_n',num2str(n)]);
else
    save(['FullSample_nlags',num2str(numlags),'_n',num2str(n)]);
end



% %========================================================================
% % Visualize fitted estimated components for the first asset
% %========================================================================

if forecast==0
    figure()
    for asset=1:n
        subplot(round(n/2)+1,round(n/2),asset)
    plot([date date(end)],reshape(Forecast.(sprintf(['TVLR_',num2str(h),'step'])).forecast_Ht(asset,asset,:),1,insample+h))
    hold on
    plot([date date(end)],reshape(Forecast.(sprintf(['TVLR_',num2str(h),'step'])).forecast_St(asset,asset,:),1,insample+h),'r')
    datetick('x','mmmyy','keeplimits','keepticks')
    legend('Total fitted volatility','Long term component')
    set(gca,'TickDir','out')
    title(['Fitted estimated volatility of: ', tickers(asset)]);
    end
end


% %========================================================================
% % Visualize results of h-step ahead predictions for assets 1 and 2
% %========================================================================
if forecast==1
    i=1;
    j=2;
    figure()
    for ii=1:length(horizons)
        h=horizons(h);
        DCC=Forecast.(sprintf(['TVLR_',num2str(h),'step'])).forecast_Ht(:,:,insample+h:end-h+1);
        CC=Vt(:,:,insample+h:end);
        dateOOS=date(:,insample+h:end);
        out=outsample+1-h;
        subplot(length(horizons),1,ii)
        plot(dateOOS,reshape(CC(i,j,:),out,1))
        hold on
        plot(dateOOS,reshape(DCC(i,j,:),out,1),'r')
        datetick('x','mmmyy','keeplimits','keepticks')
        legend('Vt','predicted Ht')
            set(gca,'TickDir','out')
        if i==j; ass=char(tickers(i)); title([num2str(h),'-step ahead Predictions: ',num2str(ass)]);
        else ass1=char(tickers(i)); ass2=char(tickers(j));
            title([num2str(h),'-step ahead Predictions: ',num2str(ass1),'-',num2str(ass2)]);
        end
    end
end
















