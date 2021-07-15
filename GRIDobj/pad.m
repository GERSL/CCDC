function DEM = pad(DEM,varargin)

%PAD add or remove a border of pixels around a GRIDobj
%
% Syntax
%
%     DEMp = pad(DEM);
%     DEMp = pad(DEM,px,val)
%
% Description
%
%     This function adds or removes a border of constant values along all 
%     edges of a GRIDobj. By default, pad adds a one-pixel wide border with
%     zeros. To control the value and amount of padding, use the arguments
%     val and px, respectively. px can be negative. In this case, the DEM
%     is cropped by the number of pixels at each grid border.
%
% Input arguments
%
%     DEM     GRIDobj
%     px      amount of padding on each of the four edges of the grid. 
%             For px=1 the resulting size of DEMp is DEM.size+2. 
%     val     pad value (default=0). This value does not have any effects 
%             if px<0.
%
% Output arguments
%
%     DEMp    padded or cropped GRIDobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEMp20 = pad(DEM,20);
%     DEMm100 = pad(DEM,-100);
%     subplot(1,2,1)
%     imageschs(DEMp20)
%     subplot(1,2,2)
%     imageschs(DEMm100)
% 
% See also: GRIDobj/crop, padarray
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. August, 2017


% check input arguments
narginchk(1,3)
p = inputParser;
p.FunctionName = 'GRIDobj/pad';

defaultval = 0;
defaultpx  = 1;
addRequired(p,'DEM')
addOptional(p,'px',defaultpx, @(x) isnumeric(x) && isscalar(x));
addOptional(p,'val',defaultval,@(x) isnumeric(x) && isscalar(x));
parse(p,DEM,varargin{:});

px  = sign(p.Results.px)*ceil(abs(p.Results.px));
val = p.Results.val;

% new size of output grid
newsize = DEM.size + px*2;

% check if size is negative and issue an error
if any(newsize < 0)
    error('the number of pixels removed from the grid edges is too large')
end

% create new array
if val == 0
    Znew = zeros(newsize,class(DEM.Z));
elseif isnan(val)
    Znew = nan(newsize,class(DEM.Z),class(DEM.Z));
else
    Znew = zeros(newsize,class(DEM.Z))+val;
end

% transfer values to new array
if px>0;
    Znew(px+1:end-px,px+1:end-px) = DEM.Z;
elseif px<0;
    abspx = abs(px);
    Znew = DEM.Z(abspx+1:end-abspx,abspx+1:end-abspx);
end

% adjust referencing
DEM.Z           = Znew;
DEM.refmat(3,:) = DEM.refmat(3,:)-px*[DEM.refmat(2,1) DEM.refmat(1,2)];
DEM.size        = size(DEM.Z);

if ~isempty(DEM.georef);
DEM.georef.SpatialRef = refmatToMapRasterReference(DEM.refmat,size(Znew));
end