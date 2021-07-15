function P = polygon2GRIDobj(DEM,MS,field)

%POLYGON2GRIDobj convert polygon to a grid
%
% Syntax
%
%     P = polygon2GRIDobj(DEM,MS)
%     P = polygon2GRIDobj(DEM,MS,pn,pv,...)
%
% Description
%
%     line2GRIDobj grids a polyline defined by a set of x and y
%     coordinates. Note that points must lie inside the grid. Segments with
%     one or two points outside the grid will not be drawn.
%     
% Input arguments
%
%     DEM    grid 
%     MS     mapstruct of a polyline as imported by shaperead. Must have
%            the fields X and Y.
%     field  field name of the numeric data to be mapped. If not
%            provided, polygon2GRIDobj will map logical values.
%
% Output arguments
%
%     P      grid (GRIDobj). Grid has the same extent and cellsize
%            as DEM.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     D = drainagebasins(FD);
%     MS = GRIDobj2polygon(D);
%     P = polygon2GRIDobj(D,MS,'ID');
%
% Note: In above example P will have a row of zeros on the top and column
% of zeros on the left side of the grid. I have not yet resolved this
% issue.
%
% See also: GRIDobj/coord2ind, GRIDobj/sub2coord, GRIDobj/getcoordinates
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. September, 2016


if nargin == 2;
    P = GRIDobj(DEM,'logical');
    writelogical = true;
    val = true;
    writeclass = @logical;
else
    P = GRIDobj(DEM,'single');
    writelogical = false;
    writeclass = @single;
end

% get DEM coordinates
[X,Y] = getcoordinates(DEM);
siz   = DEM.size;

% loop through features of mapping structure MS
for r = 1:numel(MS);
    % get coordinates
    x = MS(r).X;
    y = MS(r).Y;
    I = isnan(x) | isnan(y);
    x(I) = [];
    y(I) = [];
    
    % if the value of that attribute should be written to the grid, this
    % value is extracted here.
    if ~writelogical
        val = single([MS(r).(field)]);
    end
    
    % convert coordinates to rows and columns
    [row,col] = coord2pixel(x,y,X,Y);
    % and extract a subset of the image and apply poly2mask to that subset.
    % That is much faster than calling poly2mask for the whole image
    [BW,ext]    = getmask(row,col,siz);
    % Then, write that data back to the main grid P
    FillMat     = P.Z(ext(1):ext(2),ext(3):ext(4)); 
    FillMat(BW) = writeclass(val);
    
    P.Z(ext(1):ext(2),ext(3):ext(4)) = FillMat;

end
end

function [BW,ext] = getmask(r,c,siz)
ext = getextent(r,c);
ext = [max(floor(ext(1)),1) min(ceil(ext(2)),siz(1)) ...
       max(floor(ext(3)),1) min(ceil(ext(4)),siz(2))];
sizcrop = [ext(2)-ext(1)+1   ext(4)-ext(3)+1];
rcrop   = r - ext(1) + 1;
rcrop   = max(rcrop,1);
rcrop   = min(rcrop,sizcrop(1));

ccrop   = c - ext(3) + 1;
ccrop   = max(ccrop,1);
ccrop   = min(ccrop,sizcrop(2));

BW = poly2mask(ccrop,rcrop,sizcrop(1),sizcrop(2));
end
        
        
       

function [r,c] = coord2pixel(x,y,X,Y)

% force column vectors
x = x(:);
y = y(:);

dx  = X(2)-X(1);
dy  = Y(2)-Y(1);

c = (x-X(1))./dx + 1;
r = (y-Y(1))./dy + 1;
end


function ext = getextent(x,y)

ext = [min(x) max(x) min(y) max(y)];
end





