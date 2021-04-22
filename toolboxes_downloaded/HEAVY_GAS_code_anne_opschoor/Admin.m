
classdef Admin % in Admin.m 
methods (Static = true)

    
    %% 
    function [uLT] = UnitLT(k)
        % vult lower trianguler met enen.
        uLT= zeros(k,k);
        
        for j=1:k
            for i=j+1:k
                uLT(i,j) =1
            end
        end
         
    end
    
    
    %%
    function [Cov] = CovOneFactor(k, sigma, rho)
       % bouwt cov matrix met 1 sigma en equicorrelation
       Cov = zeros(k, k); 
       Vol = zeros(k, k); 
       Corr= eye(k, k); 
       
       for i=1:k
           Vol(i,i)= sigma;
       end
       
       
       for i=1:k
           for j=(i+1):k
              Corr(i,j)= rho;
              Corr(j,i)= Corr(i,j);
           end
       end
       
       Cov= Vol*Corr*Vol;     
    end    
    
    %%
    function [ A ] = Vech2Sym(k, vechA )
    % this function puts vech of A into A

    A= zeros(k,k);
    l=1;
        for j=1:k
            for i=j:k
                A(i,j)= vechA(l);
                A(j,i)= A(i,j);
                l=l+1;
            end
        end

    end
     
    function [ A ] = Vec2Sym(k, vecA )
    % this function puts vech of A into A

    A= zeros(k,k);
    l=1;
        for j=1:k
            for i=1:k
                A(i,j)= vecA(l);              
                l=l+1;
            end
        end

     end
    
    

    %%
    function [ vechA ] = Sym2Vech(k, A )
    % this function puts A into vech of A

    vechA= zeros(k*(k+1)/2, 1);
    l=1;
        for j=1:k
            for i=j:k
                vechA(l) = A(i,j);
                l=l+1;
            end
        end

    end

    %%
    function [ vechltA ] = LowerTr2Vech(k, ltA )
    % this function puts lower triangular A into  vech lt A

    vechltA= zeros(k*(k+1)/2, 1);
    l=1;
        for j=1:k
            for i=j:k
                vechltA(l)= ltA(i,j);
                l=l+1;
            end
        end

    end

    %%
    function [ ltA ] = Vech2LowerTr(k, vechltA )
    % this function puts vector of vech of lt A into lt A

    ltA= zeros(k,k);
    l=1;
        for j=1:k
            for i=j:k
                ltA(i,j)= vechltA(l);
                l=l+1;
            end
        end

    end
    
    %%
    function [K] = commutation(m,  n)
	
	I= eye(m*n, m*n);
	K= zeros(m*n, m*n);
	
        for j=1:n
          for i=1:m
             K((i-1)*n + j,:)= I((j-1)*m + i, :);
          end
        end
        
    end

    
    %%
    function [ D ] = duplication(k)
    % creates duplication matrix
	
	I= eye(k*(k+1)/2);
	D= zeros( k^2, k*(k+1)/2 );
	
        for j=1:k*(k+1)/2
          for i=1:k
            if (i>=j)
             D((j-1)*k + i ,:)= I((j-1)*(k - j/2) +i ,:);
             D((i-1)*k + j ,:)= I((j-1)*(k - j/2) +i ,:);
            end 
          end
        end
    
    end

	%%
    function [ L ] = elimination(k)
    % creates elimination matrix
    
    D= Admin.duplication(k); 
    L= inv(D' * D) * D'; 
	
    end
    
    
    %%
    function [DL] = duplication_gen(k)
    
    D = Admin.duplication(k);
    LT= ones(k, k);

        for i=1:k
         for j=i+1:k
           LT(i,j)= 0; 	    
         end
        end 

    LT= reshape(LT, k^2, 1) * ones(1, k*(k+1)/2 );
    DL= LT .* D;
      
    end

    %%     
    function [LL] = elimination_gen(k)

	DL= Admin.duplication_gen(k);
	LL= DL';
	
    end
    
    
    %%
    function [sqroot_matA] = sqroot_mat(A)

	[V, D]= eig(A);
    sqroot_matA= V * D.^(0.5) * inv(V);
	
    end
    
    %%
    function [sqroot_matA] = sqroot_symmat(A)

	[V, D]= eig(A);   
    sqroot_matA= V * D.^(0.5) * V';
	
    end
    
    %% 
    function [sqroot_matA] = sqroot_invsymmat(A)

	[V, D]= eig(A); %TODO: eig robust for eigen > 0
    sqroot_matA= V* (diag(diag(D).^(-1/2))) * V';
	
    end
end % end methods
end % end class