function [parameters, stderrors, LLF, ht, sk, ku, resids, summary] = garchsk(data, p, q, startingvalues, options)
%{
-----------------------------------------------------------------------
 PURPOSE:
 Leon, A., Rubio, G., and Serna, G., (2004) "Autoregressive Conditional
 Volatility, Skewness and Kurtosis"
 Brio, and Níguez, and Perote (2007) Multivariate Gram-Charlier Densities
-----------------------------------------------------------------------
 USAGE:
 [parameters, stderrors, LLF, ht, sk, ku, resids, summary] = 
 garch(data, model, distr, p, q, startingvalues, options)
 
 INPUTS:
 data:     (T x 1) vector of data
 model:    'GARCH'
 p:        positive scalar integer representing the order of ARCH
 q:        positive scalar integer representing the order of GARCH

 startingvalues: The parameter for the conditional mean is the first value
 startingvalues(1), followed by GARCH-Variance-Skewnes and Kurtosis.
 
 options:  a set of options used in fmincon
 
 OUTPUTS:
 parameters:   a vector of parameters a1, b0, b1, b2, ...
 stderrors:    standard errors estimated by the inverse Hessian (fmincon)
 LFF:          the value of the Log-Likelihood Function
 ht:           vector of conditional variance
 st:           vector of conditional skewness
 ku:           vector of conditional kurtosis
 resids:       vector of residuals
 summary:      summary of results
-----------------------------------------------------------------------
 Author:
 Alexandros Gabrielsen, alexgabriel@live.com
 Date:     05/2016
-----------------------------------------------------------------------
 Notes:
 1. Supported Specifications
    ARMAX-GARCH-SK with Time-Varying Higher Moments based on the 
    Gram-Charlier Expansion Series
 2. Some combinations of specifications and distributions may not converge.
-----------------------------------------------------------------------
%}
if nargin == 0 
    error('Data, GARCH, ARCH,  Initial Values and Options') 
end

if size(data,2) > 1
   error('Data vector should be a column vector')
end

if (length(p) > 1) | (length(q) > 1) | p < 0 | q < 0
    error('P and Q should be positive scalars')
end

% garchtype 1 = GARCH 
garchtype = 1;

% errortype: 6 = Time Varying Gram-Charlier
errortype = 6;

T=size(data,1);
m  =  max(p,q);   
% Specify Initial Values
if nargin < 6 | isempty(startingvalues)
    a1 = mean(data);
    b1 = 0.15*ones(p,1)/p;
    b2 = 0.75*ones(q,1)/q;
    b0 = (1-(sum(b1)+sum(b2)))*var(data); % unconditional
    c1 = 0.05*ones(p,1)/p;
    c2 = 0.75*ones(q,1)/q;
    c0=skewness(data);
    c1 = 0.05*ones(p,1)/p;
    c2 = 0.75*ones(q,1)/q;
    c0=(1-(sum(b1)+sum(b2)))*skewness(data); 
    d1 = 0.05*ones(p,1)/p;
    d2 = 0.75*ones(q,1)/q;
    d0=(1-(sum(b1)+sum(b2)))*(kurtosis(data)-3);  
    startingvalues = [a1; b0; b1; b2; c0; c1; c2; d0; d1; d2];    
 end 

% Specify constraints used by fmincon (A, b) and lower and upper bounds of parameters
% Example: GARCH(1,1): b0>0, b1>0, b2>0, and b1 + b2 <1, c0>0, c1>0, c2>0,
% c1 c2 < 1, d0>0, d1>0, d2>0, d1+d2<1
[r c] = size(startingvalues);
k=0;
A = [zeros(r-1,1) -eye(r-1); zeros(1,r-p-q-6) ones(1,p+q) zeros(1,6); zeros(1,r-p-q-3) ones(1,p+q) zeros(1,3); zeros(1,r-p-q) ones(1,p+q)]; 
b = [zeros(1,r-1) [1 - 1e-6] [1 - 1e-6] [1 - 1e-6]];
lowerbounds  = [-1; -1e-8*ones(size(startingvalues,1)-1-k,1)];
upperbounds  = [ 1; ones(size(startingvalues,1)-1-k,1)];
 
% define options
if nargin < 7 | isempty(options)
    options  =  optimset('fmincon');
    options  =  optimset(options , 'Algorithm ','interior-point');
    options  =  optimset(options , 'TolFun'      , 1e-006);
    options  =  optimset(options , 'TolX'        , 1e-006);
    options  =  optimset(options , 'TolCon'      , 1e-006);
    options  =  optimset(options , 'Display'     , 'iter');
    options  =  optimset(options , 'Diagnostics' , 'on');
    options  =  optimset(options , 'LargeScale'  , 'off');
    options  =  optimset(options , 'MaxIter'     , 1500);
    options  =  optimset(options , 'Jacobian'     ,'off');
    options  =  optimset(options , 'MeritFunction'     ,'multiobj');
    options  =  optimset(options , 'MaxFunEvals' , 3000);
end

[parameters, LLF, EXITFLAG, OUTPUT, Lambda, GRAD, HESSIAN] =  fmincon('garchsklik', startingvalues ,A, b, [],[],lowerbounds, upperbounds, [], options, data, garchtype, errortype, p, q, m, T);
parameters

if EXITFLAG<=0
   EXITFLAG
   fprintf(1,'Convergence has been not successful!\n')
end

[LLF,likelihoods,ht, sk, ku, resids] = garchsklik(parameters, data, garchtype, errortype, p, q, m, T);
LLF = -LLF;

% asymptotic standard errors
stderrors = sqrt(diag((HESSIAN)^(-1))); 

% t-statistics
tstats = parameters./stderrors;

% robust standard errors
h=parameters*eps;
hplus=parameters+h;
hminus=parameters-h;
likelihoodsplus=zeros(T-m,length(parameters));
likelihoodsminus=zeros(T-m,length(parameters));
for i=1:length(parameters)
    hparameters=parameters;
    hparameters(i)=hplus(i);
    [holder1,indivlike] =  garchsklik(parameters, data, garchtype, errortype, p, q, m, T);
    likelihoodsplus(:,i)=indivlike;
end
for i=1:length(parameters)
    hparameters=parameters;
    hparameters(i)=hminus(i);
     [holder1, indivlike] =  garchsklik(parameters, data, garchtype, errortype, p, q, m, T);
    likelihoodsminus(:,i)=indivlike;
end
scores=(likelihoodsplus-likelihoodsminus)./(2*repmat(h',T-m,1));
scores=scores-repmat(mean(scores),T-m,1);
S=scores'*scores;
robustSE = diag((HESSIAN^(-1))*S*(HESSIAN^(-1)));
robusttstats = parameters./robustSE;

% Saving and Organizing Results
summary.Iterations = OUTPUT.iterations;

% Statistics to be saved
A = char('Coeff', 'SErrors', 'Tstats', 'RobustSErrors', 'RobustTstats');

% Result vectors
C = char('parameters', 'stderrors', 'tstats', 'robustSE', 'robusttstats');

% Parameters to be saved   
B = strcat(char('C',...
'K',...
char(regexprep(regexp(sprintf('ARCH%1.0f/',0:p), '/', 'split'), '^.*0', '')),...
char(regexprep(regexp(sprintf('GARCH%1.0f/',0:q), '/', 'split'), '^.*0', '')),...
'', ...
'SK',...
char(regexprep(regexp(sprintf('ARCHSK%1.0f/',0:p), '/', 'split'), '^.*0', '')),...
char(regexprep(regexp(sprintf('GARCHSK%1.0f/',0:q), '/', 'split'), '^.*0', '')),...
'KU',...
char(regexprep(regexp(sprintf('ARCHKU%1.0f/',0:p), '/', 'split'), '^.*0', '')),...
char(regexprep(regexp(sprintf('GARCHKU%1.0f/',0:q), '/', 'split'), '^.*0', ''))));
   
for i = 1:size(A,1);
    for j = 1:size(B,1);
        eval(['summary.',strcat(A(i,:)),'.',strcat(B(j,:)),' =',strcat(C(i,:)),'(j);']);
    end
    clear j
end
clear i

summary.Scores=scores;
summary.LLF = LLF;
summary.AIC = -2*LLF+2*size(parameters,1);
summary.BIC = -2*LLF+size(parameters,1)*log(size(data,1));
fprintf('-------------------------------------------------\n')
fprintf('Convergence achieved after %1.0f iterations\n', OUTPUT.iterations)
fprintf('-------------------------------------------------\n')
fprintf('Parameters  Coefficients  Std Errors    T-stats\n')
fprintf('-------------------------------------------------\n')
for i = 1:size(parameters,1)
    if parameters(i) < 0
        fprintf(strcat('  %s     %1.',num2str(round(5-length(sprintf('%1.0f', abs(parameters(i)))))),'f      %1.',num2str(round(5-length(sprintf('%1.0f', stderrors(i))))),'f       %1.',num2str(round(5-length(sprintf('%1.0f',abs(tstats(i)))))),'f\n'), B(i,:), parameters(i), stderrors(i), tstats(i))
    else    
        fprintf(strcat('  %s      %1.',num2str(round(5-length(sprintf('%1.0f', abs(parameters(i)))))),'f      %1.',num2str(round(5-length(sprintf('%1.0f', stderrors(i))))),'f        %1.',num2str(round(5-length(sprintf('%1.0f', abs(tstats(i)))))),'f\n') , B(i,:), parameters(i), stderrors(i), tstats(i))
    end
end
fprintf('-------------------------------------------------\n')
fprintf('Log Likelihood: %1.0f\n', LLF)
fprintf('Akaike Information Criteron: %1.0f\n', summary.AIC)
fprintf('Bayesian Information Criteron: %1.0f\n', summary.BIC)
fprintf('-------------------------------------------------\n')

end

