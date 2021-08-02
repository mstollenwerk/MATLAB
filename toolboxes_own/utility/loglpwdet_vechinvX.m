function y = loglpwdet_vechinvX(vechinvX,powers)
invX = ivech(vechinvX);
X = inv(invX);
y = loglpwdet(X,powers);
end