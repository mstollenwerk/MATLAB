function Plot_time_series_matrix(V_mat,P_mat);

T = size(V_mat,3);
k = size(V_mat,1);
l=1;

for i=1:k
    for j=1:k
      if (i==j)  
      subplot(k,k,l), plot(squeeze(V_mat(i,j,:)));  
      elseif j<i 
      subplot(k,k,l), plot(squeeze(V_mat(i,j,:)));     
      else
       subplot(k,k,l), plot(squeeze(P_mat(i,j,:)));     
      end
      l=l+1;
    end
end