function [MS,x,y] = GRIDobj2polygon(DB,varargin)

% Conversion from drainage basin grid to polygon or polyline
%
% Syntax
%
%     MS = GRIDobj2polygon(DB)
%     MS = GRIDobj2polygon(DB,pn,pv,...)
%     [MS,x,y] = ...
%
% Description
%
%     GRIDobj2polygon converts a GRIDobj (label grid) to a mapstruct
%     that contains individual basins (or regions) as polygon or polyline
%     features. The polygone outlines or polylines run along pixel edges
%     and not pixel centers which differs from other Matlab functions to 
%     derive outlines (e.g. bwboundaries, bwperim, bwtraceboundary).
%
%     DB should be a grid of integers. Regions are those that have values
%     unequal to zero. If DB is a floating grid (single or double), then
%     GRIDobj2polygon finds unique values but discards NaNs and zero
%     values. In this case, GRIDobj2polygon returns MS with an additional
%     field 'gridval' that indicates the original value of the grid in each
%     region.
%
% Input arguments
%
%     DB    GRIDobj with drainage basins 
%     
% Parameter name/value pairs
%
%     simplify      true or {false}. Douglas-Peuker simplification. Note
%                   that for tol>0 (see tol parameter) adjacent polygons
%                   may overlap. Note that using no simplification may
%                   result in many nodes for each polygon. Use tolerance of
%                   zero if you wish that the shape is unchanged but
%                   redundant nodes are dismissed. 
%     tol           tolerance for line simplification. Zero tolerance
%                   reduces the number of nodes but leaves the geometry 
%                   unchanged. Values between 0.5 and 2 will lead to
%                   simplier shapes but usually destroy the spatial 
%                   topology of neighboring polygons.
%     minarea       scalar indicating the minimum number of pixels
%                   for a drainage basins to be included in the output.
%     geometry      'Polygon' or 'PolyLine' (or 'Line')
%     multipart     true or {false}. Enables multipart features.
%     holes         true or {false}. Enables holes. Multipart must be set
%                   to true to enable features with holes. Note that holes
%                   are currently not supported well by the algorithm.
%     waitbar       true or {false}. True will show a waitbar.
%
% Output arguments
%
%     MS            mapstruct with polygons or polylines
%     x,y           nan-separated list of coordinate vectors
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     D = drainagebasins(FD);
%     [MS,x,y] = GRIDobj2polygon(D,'Geometry','Line','simplify',true,'tol',0);
%     imageschs(DEM,DEM,'colormap','landcolor')
%     hold on
%     plot(x,y,'-r')
%
% See also: bwboundaries, bwtraceboundary, regionprops
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 13. October, 2016



% check input arguments
narginchk(1,inf)
p = inputParser;
p.FunctionName = 'GRIDobj/GRIDobj2polygon';
validgeoms  = {'Polygon','PolyLine','Line'};
addRequired(p,'DB');
addParamValue(p,'simplify',false, @(x) isscalar(x));
addParamValue(p,'tol',0, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParamValue(p,'minarea',0, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParamValue(p,'geometry','Polygon',@(x) ischar(validatestring(x,validgeoms)));
addParamValue(p,'multipart',false,@(x) isscalar(x));
addParamValue(p,'holes',false,@(x) isscalar(x));
addParamValue(p,'waitbar',false,@(x) isscalar(x));
parse(p,DB,varargin{:});

% read parsed arguments
simp  = logical(p.Results.simplify);
tol   = p.Results.tol;
minarea = p.Results.minarea;
geom  = validatestring(p.Results.geometry,validgeoms);
mp    = p.Results.multipart;
holes = p.Results.holes;
waitb = p.Results.waitbar;

% check underlying class of the grid
if isfloat(DB.Z);
    writevalue = true;
    DB2 = GRIDobj(DB,'uint32');
    I   = ~(isnan(DB.Z) | DB.Z == 0);
    [uniquevals,~,DB2.Z(I)] = unique(DB.Z(I));
    DB  = DB2;
else
    writevalue = false;
end

% identify regions and number of regions
STATS = regionprops(uint32(DB.Z),'Area','PixelIdxList');
ndb = numel(STATS);  

% go through all regions
if waitb; h = waitbar(0,'please wait'); end
    
counter = 0;
for r = 1:ndb;
    % show waitbar
    if waitb; waitbar(r/ndb,h); end
    
    % check if Area is larger than minimum area
    if STATS(r).Area <= minarea;
        continue
    else
        counter = counter+1;
    end
    
    % get subscripts of basin
    [row,col] = ind2sub(DB.size,STATS(r).PixelIdxList);
    % bw2poly returns the coordinates of the boundary
    C   = bw2poly([row col],mp,holes);
    
    % simplify line if wanted
    if simp
        if tol>0
            C = cellfun(@(x) dpsimplify(x,tol),C,'UniformOutput',false);
        else
            C = cellfun(@(x) simplify(x),C,'UniformOutput',false);
        end
    end
    
    % add nans at the end of coordinate vectors
    if numel(C)>1
        C = cellfun(@(x) [x;[nan nan]],C,'UniformOutput',false);
    end
    C = cell2mat(C);
    
    % write data to mapstruct
    MS(counter).Geometry = geom;
    [x,y] = sub2coord(DB,C(:,1),C(:,2));
    MS(counter).X = x;
    MS(counter).Y = y;
    MS(counter).ID = double(DB.Z(STATS(r).PixelIdxList(1)));
    
    if writevalue
        MS(counter).gridval = double(uniquevals(r));
    end
    
end

% create coordinate vectors if more than one output
if nargout > 1;
    for r=1:numel(MS);
        if ~isnan(MS(r).X(end))
            MS(r).X(end+1) = nan;
            MS(r).Y(end+1) = nan;
        end
    end
    
    x = {MS.X}';
    y = {MS.Y}';
    x = cell2mat(x);
    y = cell2mat(y);
end

% close waitbar
if waitb; close(h); end

end


function C = bw2poly(BW,mp,holes)

if islogical(BW);
    [r,c] = find(BW);
else
    r = BW(:,1);
    c = BW(:,2);
end
    
rc = bsxfun(@plus,r,[0 -.5 0 .5 .5 .5 0 -.5 -.5]);
cc = bsxfun(@plus,c,[0 -.5 -.5 -.5 0 .5 .5 .5 0 ]);
rc = rc(:);
cc = cc(:);

rc = rc*2;
cc = cc*2;


minrc = min(rc)-1;
maxrc = max(rc)-1;
if mp
    mincc = min(cc)-1;
else
    [mincc,ix] = min(cc);
    mincc = mincc-1;
end
maxcc = max(cc)-1;

rc = rc-minrc;
cc = cc-mincc;

siz = [maxrc-minrc+1 maxcc-mincc+1];

IX = sub2ind(siz,rc,cc);
B  = false(siz);
B(IX) = true;

if mp 
    if ~holes
        C = bwboundaries(B,4,'noholes');
        C = cellfun(@(x) modifynodelist(x),C,'UniformOutput',false);
    else
        [C,~,N] = bwboundaries(B,4,'holes');
        for iter2=1:numel(C);
            if iter2<=N
                C{iter2} = modifynodelist(C{iter2});
            else
                C{iter2} = bsxfun(@plus,C{iter2},[minrc mincc]);
                C{iter2} = C{iter2}/2;
            end
        end
            
    end
    C = C(:);
else
    C = bwtraceboundary(B,[rc(ix) cc(ix)],'S',4,inf,'counterclockwise');
    C = {modifynodelist(C)};
end
    


function nl = modifynodelist(nl)
% modify node list so that only nodes remain on pixel corners 
nl(any(mod(nl,2)==0,2),:) = [];
nl = bsxfun(@plus,nl,[minrc mincc]);
nl = nl/2;
end
end

function C = simplify(C)

% reduce number of polygon or polyline vertices
%
% This function simplifies a polyline by removing vertices whose removal
% doesn't affect the shape of the polyline. In order to work, the vertices
% must be arranged on a rectangular grid with equal dx,dy spacing.

dc  = C-circshift(C,1);
ddc = dc-circshift(dc,-1);
I   = any(ddc,2);

C   = C(I,:);

end

