function varargout = find(DEM)

%FIND Find indices of nonzero elements in GRIDobj
%
% Syntax
%
%     ix = find(DEM)
%     [r,c] = find(DEM)
%     [r,c,val] = find(DEM)
%
% Description
%
%     This function overloads the built-in function find and returns the
%     linear indices, or row and column indices of nonzero values in DEM. 
%     In contrast to the built-in find function, GRIDobj/find treats nans
%     as zero entries.
%
% Input arguments
%
%     DEM     GRIDobj
%
% Output arguments
%
%     ix      linear index
%     r       row index
%     c       column index
%     val     value of nonzero element in DEM.Z
%
% See also: find, GRIDobj/coord2ind, GRIDobj/ind2coord, GRIDobj/findcoord           
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 10. December, 2014

nargoutchk(0,3);

if isfloat(DEM.Z);
    DEM.Z(isnan(DEM.Z)) = 0;
end

if nargout == 1;
    varargout{1} = find(DEM.Z);
elseif nargout == 2;
    [varargout{1},varargout{2}] = find(DEM.Z);
else
    [varargout{1},varargout{2},varargout{3}] = find(DEM.Z);
end

    