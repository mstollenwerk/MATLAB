function [Y] = ivechchol(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cholx  = ivech(x,'lower');
Y = cholx*cholx';
end

