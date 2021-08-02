function [ K ] = comMat( r, m )
%COMMAT creates the r*m times r*m commutation matrix (see wikipedia).
%   It is defined by K*vec(X)=vec(X')

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.02.2017
%%
eye_r=eye(r);
eye_m=eye(m);
K=zeros(r*m,r*m);

for i=1:r
    for j=1:m
        K = K + kron(eye_r(:,i)*eye_m(j,:),eye_m(:,j)*eye_r(i,:));
    end
end

end

