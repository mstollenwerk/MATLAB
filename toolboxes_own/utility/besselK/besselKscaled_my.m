function y = besselKscaled_my(a,z,prec)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
str_y = besselKscaledc(a,z,prec);
y = str2double(str_y);
end

