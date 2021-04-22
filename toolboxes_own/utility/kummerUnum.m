function f = kummerUnum(a,b,x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   Tricomi's (confluent hypergeometric) function 
%   https://functions.wolfram.com/HypergeometricFunctions/HypergeometricU/27/01/
 
% doesn't work for example inputs a = 3 ,b = 2, x = 1, since 2-b (second
% input to kummer in term2) makes the function run forever.

error("Doesn't work properly for b input <=2. Explanation in code.")

term1 = gamma(1 - b)/gamma(a - b + 1)*kummer(a,b,x);
term2 = gamma(b - 1)/gamma(a)*x^(1 - b)*kummer(a - b + 1, 2 - b, x);

f = term1 - term2;

end

