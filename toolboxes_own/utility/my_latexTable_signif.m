function [ latex_string ] = my_latexTable_signif( table_, signif )
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
                if signif(i,j)==1
                    latex_string=strcat(latex_string,'$',table_{i,j}{1},'^*','$','&');
                elseif signif(i,j)==5
                    latex_string=strcat(latex_string,'$',table_{i,j}{1},'^**','$','&');
                elseif signif(i,j)==0
                    latex_string=strcat(latex_string,'$',table_{i,j}{1},'$','&');
                else
                    error('Significance Array has non-valid entries');
                end
            else
                if signif(i,j)==1
                    latex_string=strcat(latex_string,'$',table_{i,j}{1},'^*','$');
                elseif signif(i,j)==5
                    latex_string=strcat(latex_string,'$',table_{i,j}{1},'^**','$');
                elseif signif(i,j)==0
                    latex_string=strcat(latex_string,'$',table_{i,j}{1},'$');
                else
                    error('Significance Array has non-valid entries');
                end
            end
        elseif isnumeric(table_{i,j})
            if j~=col
                if signif(i,j)==1
                    latex_string=strcat(latex_string,'$',num2str(table_{i,j}),'^*','$','&');
                elseif signif(i,j)==5
                    latex_string=strcat(latex_string,'$',num2str(table_{i,j}),'^**','$','&');
                elseif signif(i,j)==0
                    latex_string=strcat(latex_string,'$',num2str(table_{i,j}),'&','$');
                else
                    error('Significance Array has non-valid entries');
                end
            else
                if signif(i,j)==1
                    latex_string=strcat(latex_string,'$',num2str(table_{i,j}),'^*','$');
                elseif signif(i,j)==5
                    latex_string=strcat(latex_string,'$',num2str(table_{i,j}),'^**','$');
                elseif signif(i,j)==0
                    latex_string=strcat(latex_string,'$',num2str(table_{i,j}),'$');
                else
                    error('Significance Array has non-valid entries');
                end                
            end
        end
    end
    latex_string=strcat(latex_string,'\\');
end

end

