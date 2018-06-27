function [DEMr,zone] = reproject2utm(DEM,res,varargin)

%REPROJECT2UTM Reproject DEM with WGS84 coordinate system to UTM-WGS84 
%
% Syntax
% 
%     [GRIDr,zone] = reproject2utm(GRID,res)
%     [GRIDr,zone] = reproject2utm(GRID,res,pn,pv,...)
%     GRIDr        = reproject2utm(GRID,GRID2)
%     GRIDr        = reproject2utm(GRID,GRID2,'method',method)
%
% Description
%
%     Reproject a grid (GRIDobj) with WGS84 geographic coordinates to UTM 
%     WGS84 (requires the mapping toolbox and image processing toolbox).
%
% Input arguments
%
%     GRID     raster (GRIDobj) with WGS84 geographic coordinates. 
%     res      spatial resolution in x- and y-direction (scalar)
%     GRID2    raster (GRIDobj) with projected coordinate system to which
%              GRID shall be projected. The resulting grid will be
%              perfectly spatially aligned (same cellsize, same upper left
%              egde, same size) with GRID2.
%     
% Parameter name/value pairs
%
%     zone       is automatically determined. If supplied, the value must
%                be a string, e.g., '32T'. Note that this function requires
%                the full grid zone reference that includes the uppercase
%                letter indicating the latitudinal band. 
%     method     interpolation method ('bilinear' (default), 'bicubic',  
%                or 'nearest')
%
% Output arguments
%
%     GRIDr    raster (GRIDobj) with UTM-WGS84 projected coordinates
%     zone     utm zone (string)
%
%
% See also: GRIDobj, imtransform, maketform, mfwdtran, minvtran, utmzone
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 12. January, 2017


% get latitude and longitude vectors
[lon,lat] = getcoordinates(DEM);

if ~isa(res,'GRIDobj');
    % and calculate centroid of DEM. The centroid is used to
    % get the utmzone
    lonc = sum(lon([1 end]))/2;
    latc = sum(lat([1 end]))/2;
    zone        = utmzone(latc,lonc);
    Nhemisphere = double(upper(zone(end)))>=78;
else
    zone = '';
    
end

% parse input arguments 
p = inputParser;
validmethods = {'bicubic','bilinear','nearest','linear'}; 
p.FunctionName = 'GRIDobj/reproject2UTM';
% required
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'res',@(x) (~isa(x,'GRIDobj') && isscalar(x) && x > 0) || isa(x,'GRIDobj'));
% optional
addParameter(p,'zone',zone,@(x) ischar(x));
addParameter(p,'method','bilinear',@(x) ischar(validatestring(x,validmethods)));

parse(p,DEM,res,varargin{:});

% get zone for output
zone = p.Results.zone;


% prepare mstruct (transformation structure) if only res supplied
if ~isa(res,'GRIDobj');
    mstruct       = defaultm('utm');
    mstruct.zone  = p.Results.zone;
    mstruct.geoid = wgs84Ellipsoid;
    mstruct       = defaultm(utm(mstruct));
    
    % use forward transformation of the corner locations of the DEM
    % to calculate the bounds of the reprojected DEM
    xlims = zeros(4,1);
    ylims = xlims;
    [xlims(1:2),ylims(1:2)]   = mfwdtran(mstruct,[min(lat) max(lat)],[min(lon) max(lon)]);
    [xlims(3:4),ylims(3:4)]   = mfwdtran(mstruct,[min(lat) max(lat)],[max(lon) min(lon)]);
    lims     = [min(xlims) max(xlims) min(ylims) max(ylims)]; 

else
    
    mstruct  = res.georef.mstruct;
    [x,y]    = getcoordinates(res);
    lims     = [min(x) max(x) min(y) max(y)];
end


% prepare tform for the image transform
T = maketform('custom', 2, 2, ...
    @FWDTRANS, ...
    @INVTRANS, ...
    []);

% calculate image transform
if ~isa(res,'GRIDobj');
    [Znew,xdata,ydata] = imtransform(flipud(DEM.Z),T,p.Results.method,...
        'Xdata',lims([1 2]),...
        'Ydata',lims([3 4]),...
        'Udata',lon([1 end]),'Vdata',lat([end 1])',...
        'XYScale',[res res],...
        'Fillvalues',nan...
        );
    % we have calculated the imtransform with 'ColumnsStartFrom' south. 
    % GRIDobjs use 'ColumnsStartFrom' north
    Znew = flipud(Znew);
    xnew = cumsum([xdata(1) repmat(res,1,size(Znew,2)-1)]);
    ynew = flipud(cumsum([ydata(1) repmat(res,1,size(Znew,1)-1)])');
else
    Znew = imtransform(flipud(DEM.Z),T,p.Results.method,...
        'Xdata',lims([1 2]),...
        'Ydata',lims([3 4]),...
        'Udata',lon([1 end]),'Vdata',lat([end 1])',...
        'XYScale',[res.cellsize res.cellsize],...
        'Fillvalues',nan...
        );
    Znew = flipud(Znew);
end



if ~isa(res,'GRIDobj');
% Construct GRIDobj
DEMr = GRIDobj(xnew,ynew,Znew);

% and include geospatial information
R = refmatToMapRasterReference(DEMr.refmat,DEMr.size);
             
% write GeoKeyDirectoryTag so that DEMr can be exported
% using GRIDobj2geotiff
DEMr.georef.SpatialRef = R;
GeoKeyDirectoryTag.GTModelTypeGeoKey = 1; % Projected coordinate system
GeoKeyDirectoryTag.GTRasterTypeGeoKey = 1; % RasterPixelIsArea
GeoKeyDirectoryTag.GTCitationGeoKey = ['PCS Name = WGS_84_UTM_zone_' mstruct.zone];
GeoKeyDirectoryTag.GeogCitationGeoKey = 'GCS_WGS_1984';
GeoKeyDirectoryTag.GeogAngularUnitsGeoKey = 9102; %Angular_Degree

% get ProjectedCSTypeGeoKey
if Nhemisphere
    hemisphere = '326';
else
    hemisphere = '327';
end

% http://www.remotesensing.org/geotiff/spec/geotiff6.html#6.3.3.1
%    WGS84 / UTM northern hemisphere:	326zz where zz is UTM zone number
%    WGS84 / UTM southern hemisphere:	327zz where zz is UTM zone number
% GeoKeyDirectoryTag.ProjectedCSTypeGeoKey = str2double([hemisphere zone(1:2)]);
GeoKeyDirectoryTag.ProjectedCSTypeGeoKey = str2double([hemisphere sprintf('%02d',str2double(zone(regexp(zone,'[0-9]'))))]); 
GeoKeyDirectoryTag.ProjLinearUnitsGeoKey = 9001; % Linear_Meter
                 
DEMr.georef.GeoKeyDirectoryTag = GeoKeyDirectoryTag;
DEMr.georef.mstruct = mstruct;
DEMr.name = [DEM.name ' (utm)'];

else
    DEMr = res;
    DEMr.Z = Znew;
    DEMr.name = [DEM.name ' (repr)'];
    zone = [];
end

% Transformation functions for imtransform
% (may want to check projfwd and projinv instead mfwdtran and minvtran)
    function x = FWDTRANS(u,~)
        [x,y] = mfwdtran(mstruct,u(:,2),u(:,1));
        x = [x y];        
    end

    function u = INVTRANS(x,~)
        [lati,long] = minvtran(mstruct,x(:,1),x(:,2));    
        u = [long lati];
    end
end

