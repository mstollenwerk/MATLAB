function y = hypergeom11_my(a,b,z,prec)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
str_y = hypergeom11c(a,b,z,prec);
y = str2double(str_y);
end

