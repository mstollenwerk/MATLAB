function y = wishpdf(X,S,dof)

k=size(S,1);
m=1;
for i=1:k
    m=m*gamma((dof+1-i)/2);
end
y=det(X)^((dof-k-1)/2)*exp(-1/2*trace(S\X))/(pi^(k*(k-1)/4)*2^(dof*k/2)*m*det(S)^(dof/2));
end