function [DEM, totalarea] = demarea(DEM,correctEdges)

%DEMAREA Calculate the corrected surface area of a DEM
%
% Syntax
%
%     [CA,totalarea] = demarea(DEM,correctEdges)
%
% Description
%
%     demarea takes a GRIDobj of elevation values and calculates the total
%     corrected surface area and the corrected surface area in each cell.
%     correctEdges is a boolean indicating whether to correct the surface
%     area in edge cells. The correction works by increasing the calculated
%     cell surface area by a factor of (1 + (8-n)/8), where n is the number
%     of neighbouring cells without valid elevation data. By default, the
%     edge correction is switched off. The surface area for each cell is
%     calculated as 1/8 of the sum of the the surface areas of the eight
%     triangles formed by the centres of the cell and two of its
%     neighbours. In the diagram below, each number represents a cell in
%     the elevation matrix, and the corrected surface area A(0) of cell 0
%     would be 1/8*Sum(A(0,1,2), A(0,2,3), A(0,3,4) ...)
%
%       1-----2-----3
%       |\    |    /|
%       |  \  |  /  |
%       |    \|/    |
%       8-----0-----4
%       |    /|\    |
%       |  /  |  \  |
%       |/    |    \|
%       7-----6-----5
%   
%
%
% Input arguments
%
%     DEM            digital elevation model (GRIDobj)
%     correctEdges   {true} or false. Determines whether to correct the 
%                    surface area in edge cells
%
% Output arguments
%
%     CA             GRIDobj that contains the area of each cell
%     totalarea      total area returned as scalar
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     CA = demarea(DEM);
%     imagesc(CA)
%
% References
%
%     This approach is based on Jenness, Jeff S. 2004. “Calculating
%     Landscape Surface Area from Digital Elevation Models.” Wildlife
%     Society Bulletin 32 (3): 829–839.
%
% Author: Jürnjakob Dugge 
% see original contribution here:
% http://www.mathworks.com/matlabcentral/fileexchange/42204-calculate-dem-surface-area
% Date: 01. April, 2014

if nargin == 1;
    correctEdges = true;
end

[totalarea, DEM.Z] = demareasub(DEM.Z, DEM.cellsize, correctEdges);
DEM.name = 'surface area';
end

function [ totalarea, cellareas ] = demareasub( Z, varargin )
%CALCAREA Calculate the corrected surface area of a DEM
%   [ totalarea, cellareas ] = CALCAREA(Z) takes a matrix of elevation
%   values and calculates the total corrected surface area and the
%   corrected surface area in each cell. Assumes that the horizontal and
%   vertical resolutions are equal.
% 
%   [ totalarea, cellareas ] = CALCAREA(Z, resolution) takes a matrix of
%   elevation values and the horizontal resolution (i.e. the distance
%   between cell centres), calculates the total corrected surface area and
%   the corrected surface area in each cell.
%
%   [ totalarea, cellareas ] = CALCAREA(Z, resolution, correctEdges) takes
%   a matrix of elevation values, the horizontal resolution, and a boolean
%   indicating whether to correct the surface area in edge cells. The
%   correction works by increasing the calculated cell surface area by a
%   factor of (1 + (8-n)/8), where n is the number of neighbouring cells
%   without valid elevation data. By default, the edge correction is
%   switched off.
%
%
%   The surface area for each cell is calculated as 1/8 of the sum of the
%   the surface areas of the eight triangles formed by the centres of the
%   cell and two of its neighbours. In the diagram below, each number
%   represents a cell in the elevation matrix, and the corrected surface
%   area A(0) of cell 0 would be 1/8*Sum(A(0,1,2), A(0,2,3), A(0,3,4) ...)
%
%   1-----2-----3
%   |\    |    /|
%   |  \  |  /  |
%   |    \|/    |
%   8-----0-----4
%   |    /|\    |
%   |  /  |  \  |
%   |/    |    \|
%   7-----6-----5
%   
%   This approach is based on Jenness, Jeff S. 2004. “Calculating Landscape
%   Surface Area from Digital Elevation Models.” Wildlife Society Bulletin
%   32 (3): 829–839.
%
%
%   Example, calculating the surface area of a pyramid with a base length
%   of 2 and a height of 1:
%   
%   [XX,YY]=meshgrid(linspace(-1,1,101),linspace(-1,1,101));
%   ZZ = 1-(max(abs(XX),abs(YY)));
%   % Estimated surface area
%   calcarea(ZZ,2/100)
%   %  Actual surface area
%   4*sqrt(2)

if length(varargin) >= 1
    resolution = varargin{1};
else
    resolution = 1;
end

if length(varargin) >= 2
    correctEdges = varargin{2};
else
    correctEdges = false;
end

Zsize = size(Z);

% Pad elevation matrix with NaN to handle the edges
B=padarray(Z,[1 1],NaN);

% Repeat the elevation matrix once for each neighbouring cell
B=repmat(B,[1 1 8]);

% Shift along principal axes first then along diagonal.
% [0 -1] shifts to the left, [1 1] shifts down and right, and so on.
shift1=[0 1 0 -1 1 1 -1 -1];
shift2=[-1 0 1 0 -1 1 1 -1];

for i = 1:8
    B(:,:,i)=circshift(B(:,:,i),[shift1(i), shift2(i)]);
end

% Trim off the edges that were introduced by the padding
B=B(2:end-1,2:end-1,:);

% Calculate elevation differences between central cells and neighbours
dZ=B-repmat(Z,[1 1 8]);

% Replace NaNs with 0, resulting in vector of zero length for NaN cells
%dZ(isnan(dZ)) = 0;

% Place the matrices next to each other
dZ = reshape(dZ, [Zsize(1), Zsize(2)*8]);


% For building the vectors describing the triangle edges,
% basically do the same that was done for the elevations (pad, shift, trim)
unitvector=ones([Zsize(1), Zsize(2),8])*resolution;
unitvector=padarray(unitvector,[1 1]);
for i = 1:8
    unitvector(:,:,i)=circshift(unitvector(:,:,i),[shift1(i), shift2(i)]);
end
unitvector = unitvector(2:end-1,2:end-1,:);

% Place the unit component matrices next to each other
unitvector = reshape(unitvector,[Zsize(1), Zsize(2)*8]);
unitvechoriz = unitvector(:,1:Zsize(2)*4);
unitvecdiag = unitvector(:,Zsize(2)*4+1:Zsize(2)*8);

% Triangle edges oriented along the grid direction,
% represented by vectors along the third dimension
vectorarrayA=cat(3,[unitvechoriz unitvechoriz],...
    zeros([Zsize(1), Zsize(2)*8]),...
    [dZ(:,1:Zsize(2)*4) dZ(:,1:Zsize(2)*4)]);

% The diagonal triangle edges (length sqrt(2))
vectorarrayB=cat(3,[unitvecdiag unitvecdiag(:,Zsize(2)*3+1:end) unitvecdiag(:,1:Zsize(2)*3)],...
    [unitvecdiag unitvecdiag(:,Zsize(2)*3+1:end) unitvecdiag(:,1:Zsize(2)*3)],...
    [dZ(:,Zsize(2)*4+1:end) dZ(:,Zsize(2)*7+1:end) dZ(:,Zsize(2)*4+1:Zsize(2)*7)]);
    

% Calculate the cross-products of the horizontal and diagonal vectors
X = cross( vectorarrayA, vectorarrayB,3);

% The area of each cell
triangleareas = reshape(sqrt(sum(X.^2,3)),[Zsize(1), Zsize(2), 8]);
triangleareas(isnan(triangleareas))=0;

cellareas = sum(0.125*triangleareas, 3);

if correctEdges
    correctionfactor = 8./sum(triangleareas ~= 0, 3);
    correctionfactor(isinf(correctionfactor))=0;
    cellareas = cellareas.*correctionfactor;
end

totalarea = sum(cellareas(:));
end