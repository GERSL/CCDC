function MASK = createmask(DEM,usehillshade)

%CREATEMASK create a binary mask using polygon mapping
%
% Syntax
%
%     MASK = createmask(DEM)
%     MASK = createmask(DEM,usehillshade)
%
% Description
%
%     createmask is an interactive tool to create a mask based on an
%     interactively mapped polygon.
%
% Input arguments
%
%     DEM           GRIDobj
%     usehillshade  use hillshade as background image ({false} or true)
%
% Output arguments
%
%     MASK   GRIDobj with logical mask
%
%
% See also: imroi, GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017

if nargin == 1;
    usehillshade = false;
end

figure
if usehillshade
   imageschs(DEM);
else
   imagesc(DEM);
end
title('Create polygon');
h = impoly;
MASK = DEM;
MASK.name = 'mask';
MASK.Z = createMask(h);
pos = getPosition(h);
delete(h);
hold on
plot(pos([1:end 1],1),pos([1:end 1],2));
hold off
title('done');
