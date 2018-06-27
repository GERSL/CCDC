function OUT = mtimes(varargin)

% overloaded multiplication for GRIDobj
%
% Note that there is no matrix multiplication defined for GRIDobj. Thus
% mtimes overloads element wise multiplication (times). 

funname = 'times';
narginchk(2,2)
OUT = builtincaller(funname,varargin{:});