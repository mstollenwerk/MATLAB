function [ varargout ] = dpcpartrans( varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin==3
    Lc=varargin{1};
    dc=varargin{2};
    rest_param=varargin{3};
    Sc=Lc*diag(dc)*Lc';
    varargout{1}=[vech(chol(Sc,'lower'))' rest_param];
elseif nargin==2
    n=varargin{1};
    Sc=ivech(varargin{2}(1:n*(n+1)/2),'lower');
    Sc=Sc*Sc';
    [Lc,dc]=dpceig(Sc);
    varargout{1}=Lc;
    varargout{2}=dc;
    varargout{3}=varargin{2}(n*(n+1)/2+1:end);
end

end

