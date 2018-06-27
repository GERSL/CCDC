function OUT = toposhielding(DEM,varargin)
%TOPOSHIELDING topographic shielding from cosmic rays
%
% Syntax
%     
%     S = toposhielding(DEM)
%     S = toposhielding(DEM,dphi,dtheta)
%
% Description
%
%     toposhielding calculates topographic shielding factors for cosmogenic
%     nuclide production rates following the method of Dunne et al. (1999).
%     The idea is to illuminate the topography from different points across
%     the sky and record for each azimuth (phi) the angle from the
%     horizontal (theta) when a pixel in the DEM is not shaded anymore.
%     Input arguments dphi and dtheta relate to the angular increments of
%     phi (0-360 deg) and theta (0-90 deg) for which shading is calculated. 
%     When called without specifying dphi and dtheta, default values of 5
%     degree are chosen. The approach is similar to the one published by 
%     Codilean (2008).
%     
%     Note that calculating topographic shielding for DEM's is
%     computationally expensive. Large DEM's may take a lot of time or even
%     cause memory issues. If topographic shielding is only needed for a
%     single point in space, try cropping the DEM to a smaller extent first
%     to save time, or measure the skyline in the field and use the
%     original formula from Dunne et al. (1999).
%     
%
% Input arguments
%
%     DEM       digital elevation model (GRIDobj)
%     dphi      horizontal angular increment (default=5°)
%     dtheta    vertical angular increment (default=5°)
%
% Output arguments
%
%     S         grid with shielding factors (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     S = toposhielding(DEM,10,5);
%     imagesc(S,[0 1])
%
%
% See also: castshadow
% 
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: May, 2014
%
% Requires the TopoToolbox v2:
% Download: http://csdms.colorado.edu/wiki/Model:TopoToolbox
% Article:  http://www.earth-surf-dynam.net/2/1/2014/esurf-2-1-2014.html
%
% When using this code, please cite the orginal work by Dunne et al. (1999) 
% and the Topotoolbox v2 (Schwanghart and Scherler, 2014).
% 
% References: 
% Codilean, A. T. (2008), Calculation of the cosmogenic nuclide production
% topographic shielding scaling factor for large areas using DEMs, Earth
% Surface Processes and Landforms, 31, 785–794.
% 
% Dunne, J., D. Elmore, and P. Muzikar (1999), Scaling factors for the 
% rates of production of cosmogenic nuclides for geometric shielding and 
% attenuation at depth on sloped surfaces, Geomorphology, 27, 3-11.
% 
% Schwanghart,W., and D. Scherler (2014), Short Communication: TopoToolbox
% 2 – MATLAB-based software for topographic analysis and modeling in Earth
% surface sciences, Earth Surf. Dynam., 2, 1–7, doi:10.5194/esurf-2-1-2014.


% Parse inputs
p = inputParser;
validationFcn = @(x) isnumeric(x) && isscalar(x) && (x > 0);
p.FunctionName = 'GRIDobj/toposhielding';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addOptional(p,'dphi',5,validationFcn);
addOptional(p,'dtheta',5,validationFcn);
parse(p,DEM,varargin{:});
dphi   = p.Results.dphi;
dtheta = p.Results.dtheta;


A = aspect(DEM);
A.Z = deg2rad(A.Z); % aspect in radians
G = gradient8(DEM,'rad'); % slopes in radians

OUT = uint8(zeros(DEM.size));
OUT = repmat(OUT,[1,1,length(0:dphi:359.99)]); % might be large depending on DEM

tic
fprintf(1,'Start image processing...\n')
n = length(0:dphi:359.99);
ct = 1;
for phi = 0 : dphi : 359.99
    fprintf(1,'  %2.0d / %2.0d ...',ct,n)
    p = deg2rad(phi);
    for theta = 90 : -dtheta : 0
        t = deg2rad(theta);

        T = OUT(:,:,ct);

        % Self shadows
        HS = acos(sin(t).*cos(G.Z) + cos(t).*sin(G.Z).*cos(p-A.Z));
        S1 = abs(HS)<pi/2;
        % Cast shadows
        SH = castshadow(DEM,phi,theta);
        S2 = SH.Z<=0.5; % not in shade
        S = min(S1,S2);
        
        ix = S>0;
        T(ix) = theta;
        OUT(:,:,ct) = T;
        
    end
    ct = ct+1;
    toc
end

% Now that we know theta, we can use Eq. 6 of Dunne et al. (1999) 
% [assuming m=2.3]
OUT = deg2rad(dphi)/(2*pi) .* sin(deg2rad(single(OUT))).^3.3;
S = sum(OUT,3);
S = 1-S;
OUT = DEM;
OUT.Z = S;









