function [mean_squared_fcsterror] = msfe(fcst_matrix,realization_matrix)

if size(fcst_matrix)~=size(realization_matrix)
    error('Input matrices have to have same size')
end
T=max(size(fcst_matrix));
mean_squared_fcsterror=NaN(T,1);
for t=1:T
    err=fcst_matrix(:,:,t)-realization_matrix(:,:,t);
    err=(err+err')./2;
    vech_=vech(err);
    mean_squared_fcsterror(t)=vech_'*vech_;
end

end