function [ data_sorted ] = countnsort( COLUMNVECTOR_data )
%sorts data ALONG ROWS in the following way: The numbers are sorted according to their 
%frequency in the data. the columns are: number appearances freqency_in_%

T=tabulate(COLUMNVECTOR_data);
[~,I]=sort(T(:,3));
data_sorted=flipud(T(I,:));

end

