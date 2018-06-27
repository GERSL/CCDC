function C = aggregate(DEMhighres,DEMlowres,varargin)

%AGGREGATE resampling a GRIDobj using aggregation
%
% Syntax
%
%     C = aggregate(A,B)
%     C = aggregate(A,B,aggfun)
%
% Description
%
%     This function resamples the grid A to C to match the extent and
%     resolution of grid B. B must spatially overlap with A and must have a
%     coarser resolution. By default, aggregate uses the mean to calculate
%     new values in grid C, but aggregate takes any other function that takes 
%     a vector and returns a scalar (e.g. median, std, ...).
%
% Input arguments
%
%     A         high resolution GRIDobj
%     B         low resolution GRIDobj. A and B must have the same 
%               coordinate system
%     aggfun    anonymous function that defines how to aggregate values. 
%               The default is @mean.
%
% Output arguments
%
%     C         GRIDobj with same extent and resolution as B
%
%
% See also: GRIDobj/resample, accumarray
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 15. June, 2017

if nargin == 2
    aggfun = @mean;
end

[x,y] = getcoordinates(DEMhighres);

IX    = zeros(DEMhighres.size);

[X,Y] = getcoordinates(DEMlowres);


sy = size(y);
for r = 1:numel(x)
    IX(:,r) = coord2ind(X,Y,repmat(x(r),sy),y);
end

z = accumarray(IX(:),DEMhighres.Z(:),[prod(DEMlowres.size) 1],aggfun);
C = DEMlowres;
C.Z = reshape(z,DEMlowres.size);
C.Z(C.Z == 0) = nan;


