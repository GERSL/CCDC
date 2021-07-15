function OUT = mpower(varargin)

% overloaded power for GRIDobj
%
% Note that there is no matrix power defined for GRIDobj. Thus
% mpower overloads element-by-element powers (power). 

funname = 'power';
narginchk(2,2)
OUT = builtincaller(funname,varargin{:});