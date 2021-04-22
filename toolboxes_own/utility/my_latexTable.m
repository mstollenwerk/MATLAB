function [ latex_string ] = my_latexTable( table_ )
%MY_LATEXTABLE should create only the body of a latex table based on a
%matlab table
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 14.12.2016

[row,col]=size(table_);

latex_string='';
for i=1:row
    for j=1:col
        if iscell(table_{i,j})
            if j~=col
                latex_string=strcat(latex_string,table_{i,j}{1},'&');
            else
                latex_string=strcat(latex_string,table_{i,j}{1});
            end
        elseif isnumeric(table_{i,j})
            if j~=col
                latex_string=strcat(latex_string,num2str(table_{i,j}),'&');
            else
                latex_string=strcat(latex_string,num2str(table_{i,j}));
            end
        end
    end
    latex_string=strcat(latex_string,'\\');
end

end

