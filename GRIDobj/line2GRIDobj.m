function L = line2GRIDobj(DEM,varargin)

%LINE2GRIDOBJ convert line to a grid
%
% Syntax
%
%     L = line2GRIDobj(DEM,x,y)
%     L = line2GRIDobj(DEM,MS)
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
%     x,y    coordinate vectors that define the vertices of a polyline.
%            Segments of the line can be separated by nans
%     MS     mapstruct of a polyline as imported by shaperead. Must have
%            the fields X and Y.
%
% Output arguments
%
%     L      gridded line (GRIDobj). Grid has the same extent and cellsize
%            as DEM.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     x = [381819  406059];
%     y = [3792813 3803583];
%     L = line2GRIDobj(DEM,x,y);
%     L = dilate(L,ones(5));
%     imageschs(DEM,L)
%
% See also: GRIDobj/coord2ind, GRIDobj/sub2coord, GRIDobj/getcoordinates
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. November, 2015


narginchk(1,3)

if nargin==2
    x = [varargin{1}.X];
    y = [varargin{1}.Y];
else
    x = varargin{1};
    y = varargin{2};
end

[rows,cols] = coord2sub(DEM,x,y);
siz = DEM.size;

L = GRIDobj(DEM,'logical');

subs = [];


for r = 2:numel(rows)
    if any(isnan(rows([r r-1])))
        continue
    end
    
    p1 = [rows(r-1) cols(r-1)];
    p2 = [rows(r)   cols(r)  ];
    
    subs = [subs;getline(p1,p2)]; %#ok<AGROW>
end

if isempty(subs)
    return
end

IX   = sub2ind(siz,subs(:,1),subs(:,2));
L.Z(IX) = true;
L.name = 'line2GRIDobj';
    

function subs = getline(p1,p2)

p = p1;
d = p2-p1;
N = max(abs(d));
s = d/N;

subs = zeros(N,2);
subs(1,:) = p1;

for ii=2:N+1
   p = p+s;
   subs(ii,:) = round(p);
end
end
end