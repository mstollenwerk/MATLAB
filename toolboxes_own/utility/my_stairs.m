function my_stairs(x,y,varargin)
%MY_STAIRS creates stair plot without vertical lines. 

[sx1,sx2] = size(x);

if sx1 > sx2
    x = x';
end

if size(y,1) == size(x,2)
    y = y';
end

if isdatetime(x)
    was_datetime = 1;
    x = datenum(x);
end

addX = arrayfun(@(x) [x NaN x],x(2:end),'uniformOutput',false);
x = [x(1) cell2mat(addX)];
if was_datetime
    x = datetime(x,'ConvertFrom','datenum');
end

y_new = [];
for ii = 1:size(y,1)
    addY = arrayfun(@(y) [y y y],y(ii,1:end-1),'uniformOutput',false);
    y_new = [y_new; cell2mat(addY) y(ii,end)];
end
plot(x,y_new,varargin{:});
end