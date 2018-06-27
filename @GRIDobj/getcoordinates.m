function [x,y] = getcoordinates(DEM)

%GETCOORDINATES get coordinate vectors of an instance of GRIDobj
%
% Syntax
%
%     [x,y] = getcoordinates(DEM)
%
% Input arguments
%
%     DEM    grid (class: GRIDobj)
%
% Output arguments
%
%     x      coordinate vector in x direction (row vector)
%     y      coordinate vector in y direction (column vector)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     [x,y] = getcoordinates(DEM);
%     surf(x,y,double(DEM.Z))
%     axis image; shading interp; camlight
%     
%
%
% See also: GRIDobj2mat
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. February, 2013

[x,y] = refmat2XY(DEM.refmat,DEM.size);