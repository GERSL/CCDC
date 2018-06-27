function DEM = excesstopography(DEM,varargin)

%EXCESSTOPOGRAPHY reconstruct surface with threshold-slope surface
% 
% Syntax
%
%    EXT = excesstopography(DEM)
%    EXT = excesstopography(DEM,pn,pv,...)
%
% Description
%
%    excesstopography uses grayscale morphological erosion to reconstruct
%    an idealised surface that features only slopes equal or less than the 
%    threshold slope. Subsequently, the function calculates the difference
%    between the DEM and idealized topography which is referred to excess
%    topography.
%
%    Note that this function takes a while to evaluate, especially if the
%    DEM and kernelsize are large. Be sure that the kernel is large enough
%    or set the option 'iterate' to true.
%
% Input arguments
%
%    DEM    Digital elevation model (class: GRIDobj)
%    
% Parameter name/value pairs
%
%    'maxgradient'   maximum gradient in degrees of slopes (default: 30°)
%    'unit'          deg (degrees = default) or tan (tangens). Applies to
%                    'maxgradient' and 'tol'.
%    'kernelsize'    side length in pixels used by the kernel. Must be 
%                    integer and odd (default: 11)
%    'output'        'difference' (default) or 'elevation'. Latter returns
%                    the eroded DEM without calculating the difference.
%    'iterate'       true (default) or false. A small kernel may not
%                    detect all excesstopography. Setting iterate to true
%                    repeatedly erodes the topography until maxgradient is
%                    reached
%    'tol'           only applicable if 'iterate' is set to true. Iteration
%                    terminates if max(gradient8(DEM))-maxgradient < tol.
%                    The default is 6e-4 ° (degrees). Note that setting tol
%                    to a very small value might cause that the algorithm
%                    uses all maxiter iterations.
%    'maxiter'       only applicable if 'iterate' is set to true. Iteration
%                    terminates if maxiter are reached. The default is 100.
%    'gradtype'      The numeric gradient can be calculated by various ways.
%                    Usually, TopoToolbox uses a two-point estimation to all 
%                    cell neighbors and takes the largest one (gradient8). 
%                    By default, excesstopography calculates the idealized 
%                    surface such that no cells have two-point gradients
%                    larger than 'maxgradient'. This option is 'G8'. 
%                    If required that the resultant surface features only
%                    gradients less than 'maxgradient' if these gradients 
%                    are calculated using central differencing (e.g. the 
%                    MATLAB built-in function gradient, than set this option
%                    to 'central'.
%                    
%
% Output arguments
%
%    EXT    excess topography (class: GRIDobj)
%
% Example
%
%    DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%    EXT = excesstopography(DEM,'maxgradient',30,'kernelsize',31);
%    imageschs(DEM,EXT)
%
% Reference
%
%    Blöthe, J.H., Korup, O., Schwanghart, W. (2015): Large landslides lie
%    low: Excess topography in the Himalaya-Karakoram ranges. Geology 43, 6, 
%    523-526.
%
% See also: GRIDobj/localtopography 
%     
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017



% Parse Inputs
p = inputParser;         
p.FunctionName = 'excesstopography';
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addParamValue(p,'maxgradient',30,@(x) isscalar(x) && (x >= 0));
addParamValue(p,'gradtype','G8', @(x) ischar(validatestring(x,{'G8','central'})));
addParamValue(p,'kernelsize',11, @(x) isscalar(x) && (rem(x,2)==1) && (x>=3));
addParamValue(p,'output','difference', @(x) ischar(validatestring(x,{'difference','elevation'})));
addParamValue(p,'iterate',true);
addParamValue(p,'maxiter',100);
addParamValue(p,'tol',6e-4);
addParamValue(p,'unit','deg', @(x) ischar(validatestring(x,{'deg','tan'})));
addParamValue(p,'verbose',false, @(x) isscalar(x));

parse(p,DEM,varargin{:});

maxgradient = p.Results.maxgradient;
gradtype    = validatestring(p.Results.gradtype,{'G8','central'});
kernelsize  = p.Results.kernelsize;
output      = validatestring(p.Results.output,{'difference','elevation'});
iterate     = p.Results.iterate;
tol         = p.Results.tol;
unit        = validatestring(p.Results.unit,{'deg','tan'});

cl          = class(DEM.Z);

% maxgradient is provided in degrees. Convert to radians
switch unit
    case 'deg'
        maxgradient = tand(maxgradient);
        tol = tand(tol);
end

switch gradtype
    case 'central'
        maxg = maxgradient;
        maxgradient = sqrt(maxgradient^2 / 2);
    otherwise
        maxg = maxgradient;
end

% prepare for ordfilt
% handle nans and edge effects
INAN   = isnan(DEM.Z);
m      = max(DEM);
DEM.Z(INAN) = inf;

% prepare kernel and offset
% rectangular kernel
domain = ones(kernelsize);
% offset 
offset = false(size(domain));
center = ceil(kernelsize/2);
offset(center,center) = true;
offset = bwdist(offset,'e');
% offset = max(offset(:))-offset;
offset = offset * DEM.cellsize * maxgradient;

% do the calculation
DEMcopy = DEM;

mG = getmaximumgradient(DEMcopy,gradtype,INAN);

if p.Results.verbose
    disp([datestr(clock) ' -- Maximum gradient: ' num2str(atand(mG))])
end

if mG > maxg;
    
    if iterate
        iter = 0;
        while ((mG-maxg) > tol) && (iter<= p.Results.maxiter);
            iter = iter+1;
            DEMcopy.Z = ordfilt2(double(DEMcopy.Z)-m,1,domain,double(offset),'zeros') + m;
            mG = getmaximumgradient(DEMcopy,gradtype,INAN);
            
            if p.Results.verbose
                disp([datestr(clock) ' -- Iteration ' num2str(iter,'%d') ', Maximum gradient: ' num2str(atand(mG))])
            end
            
        end
    else
        DEMcopy.Z = ordfilt2(double(DEMcopy.Z)-m,1,domain,double(offset),'zeros') + m;
    end
end



switch output
    case 'difference'
        DEM   = DEM-DEMcopy;
        DEM.Z = cast(DEM.Z,cl);
    case 'elevation'
        DEM.Z = cast(DEMcopy.Z,cl);
end
% replace infs with nan again
DEM.Z(INAN) = nan;
DEM.name = 'excess topography';
end

function mG = getmaximumgradient(DEM,gradtype,INAN)
DEM.Z(INAN) = nan;
switch gradtype
    case 'G8'
        G  = gradient8(DEM);
        mG = max(G);
    case 'central'
        [Gx,Gy] = gradient(DEM.Z,DEM.cellsize);
        G  = sqrt(Gx.^2 + Gy.^2);
        mG = max(G(:));
end

end