function y = iperm(x)
for ii = 1:length(x)
    y(ii) = find(x==ii);
end
end