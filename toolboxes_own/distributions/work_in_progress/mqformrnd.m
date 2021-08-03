    function [ mqform_rnd ] = mqformrnd( Sigma_, Psi_, A )
%MQFORMRND 
%
% USAGE:
%  [ mqform_rnd ] = mqformrnd( Sigma_, Psi_, A)
%
% INPUTS:
%   SIGMA_        - p by p by N array of parameter matrices.
%   PSI_          - df_ by df_ by N array of parameter matrices.
%   A             - df_ by df_ by N array of paramter matrices.
%
% OUTPUTS:
%   MQFORM_RND    - Random Matrix
%
% COMMENTS:
%   See also mvnrnd [Statistics and Machine Learning Toolbox R2019b]
%
% REFERENCES:
%      [1] 
% DEPENDENCIES:
%   mvnrnd [Statistics and Machine Learning Toolbox R2019b]

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 19.12.2019

%%
[p,~,N] = size(Sigma_);
[df_, ~, ~] = size(A);

%% Input Checking
% Are matrix dimensions valid?
if numel(size(Sigma_)) ~=2 && numel(size(Sigma_)) ~=3
    error('Sigma_ dimensions are not valid')
end
if numel(size(Psi_)) ~=2 && numel(size(Psi_)) ~=3
    error('Psi_ dimensions are not valid')
end

% Are matrices quadratic?
if size(Sigma_,1) ~= size(Sigma_,2)
    error('Sigma_ are not a quadratic matrices')
end
if size(Psi_,1) ~= size(Psi_,2)
    error('Psi_ are not a quadratic matrices')
end
if size(A,1) ~= size(A,2)
    error('A are not a quadratic matrices')
end

% Are respective matrices symmetric?
if Sigma_ ~= permute(Sigma_, [2,1,3])
    error('Sigma_ Matrices are not symmetric')
end
if Psi_ ~= permute(Psi_, [2,1,3])
    error('Psi_ Matrices are not symmetric')
end
if size(A,1) ~= size(A,2)
    error('A are not a quadratic matrices')
end


%% Random Draws
mqform_rnd = NaN(p,p,N);

for ii=1:N
    X = reshape(mvnrnd(zeros(p*df_,1), kron(Sigma_(:,:,ii), Psi_(:,:,ii))), p, df_);
    mqform_rnd(:,:,ii) = X*A(:,:,ii)*X';    
end
end
