function [ S, fcstS ] = bekkRec( param, data, p, q, initial, type )
%BEKKREC matrix series acc. to BEKK recursion.
%
% USAGE:
%  [S, fcstS] = bekkRec(param,data,p,q,initial,type)
%
% INPUTS:
%   PARAM    - Vector of parameters.  The form of the parameters depends on the TYPE.
%   DATA     - K by K by T array of covariance estimators or cell of length 
%              T of K by K covariance estimators (e.g. realized covariance).
%   P        - Number of lagged covariance matrix data.
%   Q        - Number of lagged recursion matrices, i.e. conditional covariance matrices.
%   INITIAL  - [Optional] Initialization matrix of the recursion. Default is time series average
%   TYPE     - [Optional] String, one of 'Scalar' (Default) ,'Diagonal' or 'Full'.
%
% OUTPUTS:
%   S        - K by K by T dimension matrix of conditional covariance matrices.
%   fcstS    - K by K array, 1-step ahead covariance matrix forecast.
%
% COMMENTS:
%   
%   PARAMETER INPUT FORMAT:
%       'Scalar':
%       [c a(1) ... a(p) b(1) ... b(q)]'  (all scalars)
%       'Diagonal'
%       [c diag(A(:,:,1))' ... diag(A(:,:,p))' diag(B(:,:,1))' ... diag(B(:,:,q))']'
%       'Full'
%       [c f(A(:,:,1)) ... f(A(:,:,p)) f(B(:,:,1)) ... f(B(:,:,q))]'
%       where c = vech(chol(C,'lower'))' and f(M) = M(:)'
%
% REFERENCES:
%      [1] Engle and Kroner (1995) - Multivariate simultaneous 
%          generalized arch.
%
% based on:
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012
%
% See also

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 13.02.2017

if isnumeric(data) % Double input  
    [k,~,T]=size(data);
elseif iscell(data) % Cell input  
    k = size(data{1},1);
    T = max(size(data));
end
[C,A,G,B] = bekk_parameter_transform(param,p,0,q,k,type);
%% Recursions either cell or 3d array (cell is usually faster)
if isnumeric(data) % Double input
    S = NaN(k,k,T+1);
    % Default Initialization
    if isempty(initial)
        initial = mean(data,3);    
    end    
    for i=1:T+1
        S(:,:,i) = C;
        for j=1:p
            if (i-j)<=0
                S(:,:,i) = S(:,:,i) + A(:,:,j)'*initial*A(:,:,j);
            else
                S(:,:,i) = S(:,:,i) + A(:,:,j)'*data(:,:,i-j)*A(:,:,j);
            end
        end   
        for j=1:q
            if (i-j)<=0
                S(:,:,i) = S(:,:,i) + B(:,:,j)'*initial*B(:,:,j);
            else
                S(:,:,i) = S(:,:,i) + B(:,:,j)'*S(:,:,i-j)*B(:,:,j);
            end
        end
    end
    fcstS = S(:,:,end);
    S = S(:,:,1:end-1);
elseif iscell(data) % Cell input  
    S = cell(1,T+1);
    % Default Initialization
    if isempty(initial)
        initial = mean(cat(3,data{:}),3);    
    end
    for i=1:T+1
        S{i} = C;
        for j=1:p
            if (i-j)<=0
                S{i} = S{i} + A(:,:,j)'*initial*A(:,:,j);
            else
                S{i} = S{i} + A(:,:,j)'*data{i-j}*A(:,:,j);
            end
        end   
        for j=1:q
            if (i-j)<=0
                S{i} = S{i} + B(:,:,j)'*initial*B(:,:,j);
            else
                S{i} = S{i} + B(:,:,j)'*S{i-j}*B(:,:,j);
            end
        end
    end
    fcstS = S{end};
    S = S{1:end-1};
else error('data must be 3-d array of cell')
end

end

