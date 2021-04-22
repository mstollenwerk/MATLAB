function latex = latex_inSample_intercept(INPUT,n,Vt,TVLRHEAVY)
% These settings are used if the corresponding optional inputs are not given.
% Pivoting of the input data swithced off per default:
% Default mode for applying input.tableDataFormat:
if ~isfield(INPUT,'dataFormatMode'),INPUT.dataFormatMode = 'column';end
% Sets the default display format of numeric values in the LaTeX table to '%.4f'
% (4 digits floating point precision).
if ~isfield(INPUT,'dataFormat'),INPUT.dataFormat = {'%.3f'};end
% Define what should happen with NaN values in input.tableData:
if ~isfield(INPUT,'nanString'),INPUT.dataNanString = '-';end
% Specify the alignment of the columns:
% 'l' for left-justified, 'c' for centered, 'r' for right-justified
if ~isfield(INPUT,'tableColumnAlignment'),INPUT.tableColumnAlignment ='c';end
% Specify whether the table has borders:
% 0 for no borders, 1 for borders
if ~isfield(INPUT,'tableBorders'),INPUT.tableBorders = 0;end
% Other optional fields:
if ~isfield(INPUT,'tableCaption'),INPUT.tableCaption = 'Unconditional covariance matrix and estimated intercept parameters';end
if ~isfield(INPUT,'tableLabel'),INPUT.tableLabel = 'int_estim';end
if ~isfield(INPUT,'makedocument'),INPUT.makedocument='1'; end
%====================================================================
% Author : Manuela Braione, CORE 2016
%====================================================================


%total number of columns and rows
numberCols= n*2;
numberRows= n;

%compute unconditional covariance matrix
A=zeros(n,n);
for t=1:size(Vt,3)
    a=Vt(:,:,t);
    A=A+a;
end
P0=A/size(Vt,3);


%%%% obtain cell array for the table data and labels
C = cell(numberRows, numberCols);
Q = cell(numberRows, numberCols);

%define length of decimal numbers
dataFormatArray = cell(1,2);
dataFormatArray(1)=INPUT.dataFormat;
dataFormatArray(2)={'%.0f'};

% make table header lines:
hLine= '\hline';
cLine = ['\cline{', num2str(ceil(numberCols/2-1+1)),'-',num2str(ceil(numberCols/2+1)),'}'];
if INPUT.tableBorders
    header = ['\begin{tabular}{|',repmat([INPUT.tableColumnAlignment,'|'],1,size(C,2)),'}'];
else
    header = ['\begin{tabular}{',repmat(INPUT.tableColumnAlignment,1,n),'@{\hskip 0.4in}',repmat(INPUT.tableColumnAlignment,1,n),'}'];
end
capt=['\caption{',INPUT.tableCaption,'}'];
beforeHeader='\resizebox{.6\textwidth}{!}{';
beginning = {'\begin{table}[h!]';capt;'\centering';'\smallskip';beforeHeader;header;hLine};
latex=[];
latex=[latex; beginning];
latex(end+1)={['\multicolumn{',num2str(n),'}{c}{\bf{E($V_t$)}} & \multicolumn{',num2str(n),'}{c}{$\bf{\bar{\Lambda}}$} \\']};
latex(end+1)={hLine};


% generate table
Lambdabar     = ivech(TVLRHEAVY.int')*ivech(TVLRHEAVY.int')';
Lambda_stderr = math(TVLRHEAVY.std(1:n*(n+1)/2)');

%%%% put in main table
for i=1:size(Vt,1)
    C(i,1:end)=num2cell([P0(i,:) Lambdabar(i,:)]);
    Q(i,1:end)=num2cell([NaN(1,n) Lambda_stderr(i,:)]);
end


for i=1:size(Vt,1)
    riga=C(i,:);
    riga2=Q(i,:);
    
    for j=1:length(riga)
        el=riga(j);
        el2=riga2(j);
        if iscellstr(el), el=char(el); end
        if iscell(el), el=cell2mat(el); end
                if iscellstr(el2), el2=char(el2); end
        if iscell(el2), el2=cell2mat(el2); end
        dataValue = num2str(el,dataFormatArray{1});
        dataValue2 = num2str(el2,dataFormatArray{1});
        if j==1; Riga=[dataValue] ; else
            if strcmp(dataValue2,'NaN');   
                Riga =[Riga,' & ', '$\displaystyle\mathop{',dataValue,'}_{' '}$'];
        else
                Riga =[Riga,' & ', '$\displaystyle\mathop{',dataValue,'}_{(', dataValue2,')}$'];
            end
        end
    end
    latex(end+1)={[Riga,' \\']};
end
latex(end+1)={hLine};



% make table footer lines:
tableLabel='Ratio_comp';

footer = {'\end{tabular}'; ...
    ['\label{tab_',tableLabel,'}'];'}';...
    '\end{table}'};
    latex = [latex;footer];

if INPUT.makedocument
    latexHeader = {'\documentclass[a4paper,10pt]{article}';'\usepackage{graphicx}';'\begin{document}'};
    latex = [latexHeader;latex];
    latexFooter = {'\end{document}'};
    latex=[latex; latexFooter];
    latex=char(latex);
else
    latex=char(latex);
end
end

function res = math(v)
% Authors: Juan-Pablo Ortega, Lyudmila Grigoryeva

% Inverse of the vech operator
% Given a (vertical) 1/2(n+1)n dimensional vector v, z is the corresponding
% nxn symmetric matrix obtained by filling z by rows

u = v';
n = round(-.5 + .5 * sqrt(1 + 8 * length(v)));

res = zeros(n);

j = 1;
for i = 1:n
    res(i, i:n) = u(j:j + n - i);
    j = j + n - i + 1;
end

res = res + res' - diag(diag(res));
end