function h = imagesc(DEM,varargin)

%IMAGESC Scale data in GRIDobj and display as image object
%
% Syntax
%
%     h = imagesc(DEM,varargin)
%
% Description
%
%     see imagesc 
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. January, 2013

[x,y] = refmat2XY(DEM.refmat,DEM.size);
ht = imagesc(x,y,DEM.Z,varargin{:});

axis xy
axis image

if nargout == 1;
    h = ht;
end
