function OUT = mrdivide(varargin)

% overloaded right division for GRIDobj
%
% Note that there is no matrix right division defined for GRIDobj. Thus
% mrdivide overloads element wise right division (rdivide). 

funname = 'rdivide';
narginchk(2,2)
OUT = builtincaller(funname,varargin{:});