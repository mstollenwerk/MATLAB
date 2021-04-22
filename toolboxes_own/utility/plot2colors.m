function plot2colors( x, y, mark_id )
%PLOT2COLORS plots x against y using red color when mark_id is 1.

%% Input Check
if any(mark_id(mark_id~=0)~=1)
    error('mark_id must be vector of zeros and ones');
end

%%
mark_id(mark_id==0) = NaN;
y_red = y.*mark_id;
plot(x,y)
hold on
plot(x,y_red,'r')

end
