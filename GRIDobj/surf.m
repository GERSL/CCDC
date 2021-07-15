function ht = surf(DEM,varargin)

%SURF surface plot for GRIDobj
%
% Syntax
%
%     surf(DEM)
%     surf(DEM,A)
%     surf(...,pn,pv,...)
%     h = ...
%
% Description
%
%     surf for GRIDobj overloads the surf command and thus provides fast
%     access to 3D visualization of digital elevation models. Note that
%     GRIDobj/surf automatically resamples the DEM so that the maximum of
%     rows or columns does not exceed 1000. This ensures that the surface
%     is efficiently drawn.
%
% Input arguments
%
%     DEM          digital elevation model (GRIDobj)
%     A            grid to define color (GRIDobj) 
%
% Parameter name/values
%
%     'exaggerate'   height exaggeration, default = 1
%      
%     and all property name/value pairs allowed by surface objects.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     surf(DEM)
%     camlight
%
% See also: imageschs
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. January, 2013


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change threshold (maximum number of rows or cols) here 
maxrowsorcols = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if max(DEM.size)>maxrowsorcols
    
    resamplecellsize = max(DEM.size)/(maxrowsorcols-1) * DEM.cellsize;
    DEM = resample(DEM,resamplecellsize);
    
    if ~isempty(varargin) && isa(varargin{1},'GRIDobj')
        A = varargin{1};
        varargin(1) = [];
        A = resample(A,DEM);
        overlay = true;
    else
        overlay = false;
    end
         
    
else
    if ~isempty(varargin) && isa(varargin{1},'GRIDobj')
        A = varargin{1};
        varargin(1) = [];
        validatealignment(DEM,A);
        overlay = true;
    else
        overlay = false;
    end
end

% go through varargin to find 'exaggerate'
TF = strcmpi('exaggerate',varargin);
ix = find(TF);
if ~isempty(ix)
    exagg = varargin{ix+1};
    varargin([ix ix+1]) = [];
else
    exagg = 1;
end

% go through varargin to find 'block'
TF = strcmpi('block',varargin);
ix = find(TF);
if ~isempty(ix)
    block = varargin{ix+1};
    varargin([ix ix+1]) = [];
else
    block = false;
end

% go through varargin to find 'baselevel'
TF = strcmpi('baselevel',varargin);
ix = find(TF);
if ~isempty(ix)
    baselevel = varargin{ix+1};
    varargin([ix ix+1]) = [];
else
    baselevel = min(DEM)*.8;
end

[x,y] = refmat2XY(DEM.refmat,DEM.size);    

if overlay
    h = surf(x,y,double(DEM.Z),double(A.Z),varargin{:});
else
    h = surf(x,y,double(DEM.Z),varargin{:});
end

exaggerate(gca,exagg);
shading interp
% camlight

if block 
    
    if any(isnan(DEM));
        error('TopoToolbox:wronginput','DEM must not have NaNs')
    end
    facecolor  = [.6 .6 .6];
    
    xp = x(:);
    xp = [xp(1); xp; xp(end:-1:1)];
    yp = repmat(y(1),size(xp));
    zp = DEM.Z(1,:);
    zp = zp(:);
    zp = [baselevel;zp;repmat(baselevel,size(zp))];
    p(1) = patch(xp,yp,zp,facecolor,'EdgeColor','none','FaceColor',facecolor);
    
    yp = repmat(y(end),size(xp));
    zp = DEM.Z(end,:);
    zp = zp(:);
    zp = [baselevel;zp;repmat(baselevel,size(zp))];
    p(2) = patch(xp,yp,zp,facecolor,'EdgeColor','none','FaceColor',facecolor);
    
    yp = y(:);
    yp = [yp(1); yp; yp(end:-1:1)];
    xp = repmat(x(1),size(yp));
    zp = DEM.Z(:,1);
    zp = zp(:);
    zp = [baselevel;zp;repmat(baselevel,size(zp))];
    p(3) = patch(xp,yp,zp,facecolor-0.1,'EdgeColor','none','FaceColor',facecolor-.1);
    
    xp = repmat(x(end),size(yp));
    zp = DEM.Z(:,end);
    zp = zp(:);
    zp = [baselevel;zp;repmat(baselevel,size(zp))];
    p(4) = patch(xp,yp,zp,facecolor-0.1,'EdgeColor','none','FaceColor',facecolor-.1);
    
    set(p,'FaceLighting','none')
    
    
    axislim = axis;
    axislim(5) = baselevel;
    axis(axislim)
    
    
    
%     
%     % plot patches on to specified depth
%     B = bwboundaries(~isnan(DEM.Z));
%     for r = 1:numel(B);
%         [xb,yb] = sub2coord(DEM,B{r}(:,1),B{r}(:,2));
%         zb = DEM.Z(sub2ind(DEM.size,B{r}(:,1),B{r}(:,2)));
%         xb = [xb(1);xb; xb(end:-1:1)];
%         yb = [yb(1);yb; yb(end:-1:1)];
%         zb = [bottomelev;zb; zeros(numel(zb),1)+bottomelev];
%         p  = patch(xb,yb,zb,zeros(size(zb)));
%         set(p,'FaceColor',[0.7 0.7 0.7]);
%     end
end
        


    
if nargout == 1
    ht = h;
end


function exaggerate(axes_handle,exagfactor)

% elevation exaggeration in a 3D surface plot
%
% Syntax
%
%     exaggerate(axes_handle,exagfactor)
%
% Description
%
%     exaggerate is a simple wrapper for calling set(gca...). It controls
%     the data aspect ratio in a 3D plot and enables elevation
%     exaggeration.   
%
% Input
%
%     axes_handle   digital elevation model
%     exagfactor    exaggeration factor (default = 1)
%
% Example
%
%     load exampledem
%     for r = 1:4;
%     subplot(2,2,r);
%     surf(X,Y,dem); exaggerate(gca,r);
%     title(['exaggeration factor = ' num2str(r)]);
%     end
%
% 
% See also: SURF
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 6. September, 2010

if nargin == 1;
    exagfactor = 1;
end
axis(axes_handle,'image');
set(axes_handle,'DataAspectRatio',[1 1 1/exagfactor]);



