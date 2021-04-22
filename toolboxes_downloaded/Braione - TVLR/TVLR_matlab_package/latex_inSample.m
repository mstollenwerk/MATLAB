function latex = latex_inSample(latex_tab,n,TVLRHEAVY)
%==========================================================================
% author : Manuela Braione, CORE 2016
%==========================================================================
% These settings are used if the corresponding optional inputs are not given.
% Sets the default display format of numeric values in the LaTeX table to '%.3f'
% (3 digits floating point precision).
if ~isfield(latex_tab,'dataFormat'),latex_tab.dataFormat = {'%0.2f'};end
if ~isfield(latex_tab,'tableCaption'),latex_tab.tableCaption = 'Full-sample estimation results';end
if ~isfield(latex_tab,'tableLabel'),latex_tab.tableLabel = 'FS_estim';end
if ~isfield(latex_tab,'makedocument'),latex_tab.makedocument='1'; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find number of assets
if n<6; Latex_table='disp_all';
    N=n*(n+1)/2;
    numberCols  = 2;
    numberGARCH = 2*n ;
    numberRows  = numberGARCH+4+3;
else Latex_table='disp_short';
    N=n*(n+1)/2;
    numberCols  = 2;
    numberGARCH = 2;
    numberRows  = numberGARCH+4+3;
end

%define columns and rows labels
rowlabel_1=[];
rowlabel_2=[];
if strcmp(Latex_table,'disp_all')
    for i=1:n
        rowlabel1=(sprintf(['gamma_',num2str(i)]));
        rr={['$\',rowlabel1,'$']};
        rowlabel2=(sprintf(['delta_',num2str(i)]));
        rr2={['$\',rowlabel2,'$']};
        rowlabel_1=[rowlabel_1 rr];
        rowlabel_2=[rowlabel_2 rr2];
    end
    rowlabel=[rowlabel_1 rowlabel_2];
    rowlabel=[rowlabel,{'$\theta$','$\omega$','$\alpha$','$\beta$'}];
    
else
    rowlabel={'$\bar \gamma$ ','$\bar \delta$','$\theta$','$\omega$','$\alpha$','$\beta$'};
    
end

rowLabel_bis={'Loglik','AIC','BIC'};

%%%% obtain cell array for the table data and labels
C = cell(numberRows, numberCols);
C(1:end-3,1)=rowlabel;
C(end-2:end,1)=rowLabel_bis;


%define length of decimal numbers
dataFormat = cell(1,2);
dataFormat(1)=latex_tab.dataFormat;
dataFormat(2)={'%0.0f'};


% make table header lines:
hLine= '\hline';
cLine = ['\cline{', num2str(ceil(numberCols/2-1+1)),'-',num2str(ceil(numberCols/2+1)),'}'];
header = ['\begin{tabular}{l',repmat('c',1,size(C,2)-1),'}'];

capt=['\caption{',latex_tab.tableCaption,'.}'];
beforeHeader='\resizebox{.8\textwidth}{!}{';
beginning = {'\begin{table}[h!]';capt;'\centering';'\smallskip';beforeHeader;header;hLine};
latex=[];
latex=[latex; beginning];
latex(end+1)={'& \bf{TVLR-HEAVY} \\'};
latex(end+1)={hLine};


% generate table
if strcmp(Latex_table,'disp_all')
    alpha      =TVLRHEAVY.DCC(1);
    beta       =TVLRHEAVY.DCC(2);
    teta       =TVLRHEAVY.Longrun(1);
    omega      =TVLRHEAVY.Longrun(2);
    gamma      =TVLRHEAVY.GARCH(:,1);
    delta      =TVLRHEAVY.GARCH(:,2);
    totLL      =TVLRHEAVY.loglik;
    aic        =TVLRHEAVY.aic_bic(1);
    bic        =TVLRHEAVY.aic_bic(2);
    
    C(1:n,2)    =num2cell(gamma);
    C(n+1:2*n,2)=num2cell(delta);
    C(end-6,2)=num2cell(teta);
    C(end-5,2)=num2cell(omega);
    C(end-4,2)=num2cell(alpha);
    C(end-3,2)=num2cell(beta);
    C(end-2,2)=num2cell(totLL);
    C(end-1,2)=num2cell(aic);
    C(end,2)=num2cell(bic);
    
    
    if isfield(TVLRHEAVY,'std')
        Q=C;
        pval_alpha =TVLRHEAVY.std(end-1);
        pval_beta  =TVLRHEAVY.std(end);
        pval_teta  =TVLRHEAVY.std(N+1);
        pval_omega =TVLRHEAVY.std(N+2);
        pval_gamma =TVLRHEAVY.std(N+2+1:N+2+n);
        pval_delta =TVLRHEAVY.std(N+2+n+1:N+2+2*n);
        
        Q(1:n,2)    =num2cell(pval_gamma);
        Q(n+1:2*n,2)=num2cell(pval_delta);
        Q(end-6,2)=num2cell(pval_teta);
        Q(end-5,2)=num2cell(pval_omega);
        Q(end-4,2)=num2cell(pval_alpha);
        Q(end-3,2)=num2cell(pval_beta);
        
        
        for i=1:size(C,1)-3
            riga=C(i,:);
            riga2=Q(i,:);
            
            for j=1:length(riga)
                el=riga(j);
                el2=riga2(j);
                if iscellstr(el), el=char(el); end
                if iscell(el), el=cell2mat(el); end
                if iscellstr(el2), el2=char(el2); end
                if iscell(el2), el2=cell2mat(el2); end
                dataValue = num2str(el,dataFormat{1});
                dataValue2 = num2str(el2,dataFormat{1});
                if j==1; Riga=['\bf{',dataValue,'}'] ; else
                    if strcmp(dataValue,'NaN'); Riga =[Riga,' & ', '-']; else
                        Riga =[Riga,' & ', '$\displaystyle\mathop{',dataValue,'}_{(', dataValue2,')}$']; end
                end
            end
            latex(end+1)={[Riga,' \\']};
        end
        latex(end+1)={['\hline']};
        
        for i=size(C,1)-2:size(C,1)
            riga=C(i,:);
            for j=1:length(riga)
                el=riga(j);
                if iscellstr(el), el=char(el); end
                if iscell(el), el=cell2mat(el); end
                dataValue = num2str(el,dataFormat{2});
                if j==1; Riga=['\bf{',dataValue,'}']; else
                    if strcmp(dataValue,'NaN'); Riga =[Riga,' & ', '-']; else
                        Riga =[Riga,' & ', dataValue]; end; end
            end
            latex(end+1)={[Riga,' \\']};
        end
        latex(end+1)={['\hline']};
        
    else
        for i=1:size(C,1)-3
            riga=C(i,:);
            
            for j=1:length(riga)
                el=riga(j);
                if iscellstr(el), el=char(el); end
                if iscell(el), el=cell2mat(el); end
                dataValue = num2str(el,dataFormat{1});
                if j==1; Riga=['\bf{',dataValue,'}'] ; else
                    if strcmp(dataValue,'NaN'); Riga =[Riga,' & ', '-']; else
                        Riga =[Riga,' & ',dataValue]; end
                end
            end
            latex(end+1)={[Riga,' \\']};
        end
        latex(end+1)={['\hline']};
        
        for i=size(C,1)-2:size(C,1)
            riga=C(i,:);
            for j=1:length(riga)
                el=riga(j);
                if iscellstr(el), el=char(el); end
                if iscell(el), el=cell2mat(el); end
                dataValue = num2str(el,dataFormat{2});
                if j==1; Riga=['\bf{',dataValue,'}']; else
                    if strcmp(dataValue,'NaN'); Riga =[Riga,' & ', '-']; else
                        Riga =[Riga,' & ', dataValue]; end; end
            end
            latex(end+1)={[Riga,' \\']};
        end
        latex(end+1)={['\hline']};
    end
    
else
    
    alpha      =TVLRHEAVY.DCC(1);
    beta       =TVLRHEAVY.DCC(2);
    teta       =TVLRHEAVY.Longrun(1);
    omega      =TVLRHEAVY.Longrun(2);
    gamma      =mean(TVLRHEAVY.GARCH(:,1));
    delta      =mean(TVLRHEAVY.GARCH(:,2));
    totLL      =TVLRHEAVY.loglik;
    aic        =TVLRHEAVY.aic_bic(1);
    bic        =TVLRHEAVY.aic_bic(2);
    
    C(1,2)    =num2cell(gamma);
    C(2,2)    =num2cell(delta);
    C(end-6,2)=num2cell(teta);
    C(end-5,2)=num2cell(omega);
    C(end-4,2)=num2cell(alpha);
    C(end-3,2)=num2cell(beta);
    C(end-2,2)=num2cell(totLL);
    C(end-1,2)=num2cell(aic);
    C(end,2)=num2cell(bic);
    
    
    if ~isfield(TVLRHEAVY,'std')
        Q=C;
        pval_alpha =TVLRHEAVY.std(end-1);
        pval_beta  =TVLRHEAVY.std(end);
        pval_teta  =TVLRHEAVY.std(N+1);
        pval_omega =TVLRHEAVY.std(N+2);
        pval_gamma =mean(TVLRHEAVY.std(N+2+1:N+2+n));
        pval_delta =mean(TVLRHEAVY.std(N+2+n+1:N+2+2*n));
        
        Q(1,2)    =num2cell(pval_gamma);
        Q(2,2)    =num2cell(pval_delta);
        Q(end-6,2)  =num2cell(pval_teta);
        Q(end-5,2)  =num2cell(pval_omega);
        Q(end-4,2)  =num2cell(pval_alpha);
        Q(end-3,2)  =num2cell(pval_beta);
        
        
        for i=1:size(C,1)-3
            riga=C(i,:);
            riga2=Q(i,:);
            
            for j=1:length(riga)
                el=riga(j);
                el2=riga2(j);
                if iscellstr(el), el=char(el); end
                if iscell(el), el=cell2mat(el); end
                if iscellstr(el2), el2=char(el2); end
                if iscell(el2), el2=cell2mat(el2); end
                dataValue = num2str(el,dataFormat{1});
                dataValue2 = num2str(el2,dataFormat{1});
                if j==1; Riga=['\bf{',dataValue,'}'] ; else
                    if strcmp(dataValue,'NaN'); Riga =[Riga,' & ', '-']; else
                        Riga =[Riga,' & ', '$\displaystyle\mathop{',dataValue,'}_{(', dataValue2,')}$']; end
                end
            end
            latex(end+1)={[Riga,' \\']};
        end
        latex(end+1)={['\hline']};
        
        for i=size(C,1)-2:size(C,1)
            riga=C(i,:);
            for j=1:length(riga)
                el=riga(j);
                if iscellstr(el), el=char(el); end
                if iscell(el), el=cell2mat(el); end
                dataValue = num2str(el,dataFormat{2});
                if j==1; Riga=['\bf{',dataValue,'}']; else
                    if strcmp(dataValue,'NaN'); Riga =[Riga,' & ', '-']; else
                        Riga =[Riga,' & ', dataValue]; end; end
            end
            latex(end+1)={[Riga,' \\']};
        end
        latex(end+1)={['\hline']};
        
    else
        for i=1:size(C,1)-3
            riga=C(i,:);
            
            for j=1:length(riga)
                el=riga(j);
                if iscellstr(el), el=char(el); end
                if iscell(el), el=cell2mat(el); end
                dataValue = num2str(el,dataFormat{1});
                if j==1; Riga=['\bf{',dataValue,'}'] ; else
                    if strcmp(dataValue,'NaN'); Riga =[Riga,' & ', '-']; else
                        Riga =[Riga,' & ',dataValue]; end
                end
            end
            latex(end+1)={[Riga,' \\']};
        end
        latex(end+1)={['\hline']};
        
        for i=size(C,1)-2:size(C,1)
            riga=C(i,:);
            for j=1:length(riga)
                el=riga(j);
                if iscellstr(el), el=char(el); end
                if iscell(el), el=cell2mat(el); end
                dataValue = num2str(el,dataFormat{2});
                if j==1; Riga=['\bf{',dataValue,'}']; else
                    if strcmp(dataValue,'NaN'); Riga =[Riga,' & ', '-']; else
                        Riga =[Riga,' & ', dataValue]; end; end
            end
            latex(end+1)={[Riga,' \\']};
        end
        latex(end+1)={['\hline']};
    end
end



% make table footer lines:
footer = {'\end{tabular}'; ...
    ['\label{tab_',latex_tab.tableLabel,'}'];'}';'\end{table}'};
    latex = [latex;footer];
if latex_tab.makedocument==1
    latexHeader = {'\documentclass[a4paper,10pt]{article}';'\usepackage{graphicx}';'\begin{document}'};
    latex = [latexHeader;latex];
    latexFooter = {'\end{document}'};
    latex=[latex; latexFooter];
end
latex=char(latex);
end