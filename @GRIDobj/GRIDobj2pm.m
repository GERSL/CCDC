function [X,ix,txt,txtname] = GRIDobj2pm(varargin)

%GRIDobj2pm combine several GRIDobj into a predictor matrix
%
% Syntax
%
%     [X,ix,txt,txtname] = GRIDobj2pm(A,B,C,...)
%
% Description
%
%     GRIDobj2pm transforms several instances of GRIDobj to a predictor
%     matrix X. The function removes rows with nans. X(:,m) = A.Z(ix) 
%     where A is the m-th GRIDobj in the input argument list and ix
%     is the linear index that refers to non-nan values in A.Z.
%
%     Note that the predictor matrix does not include an offset, i.e. a
%     column with ones.
%
% Input arguments
%
%     A,B,C,...   spatially aligned GRIDobj with numeric values
%
% Output arguments
%
%     X           predictor matrix
%     ix          index into GRIDobj
%     txt         cell array with variable names of input arguments
%     txtname     cell array with strings obtained from A.name, B.name, ...
%
% Example: A simple landform classification using kmeans
%
%     DEM  = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD   = FLOWobj(DEM,'preprocess','carve');
%     G    = gradient8(DEM);
%     logA = log(flowacc(FD));
%     C    = curvature(DEM);
%     [X,ix,txt,txtname] = GRIDobj2pm(G,C,logA);
%     X = zscore(X);
%     IDX = kmeans(X,5);
%     CL = GRIDobj(DEM);
%     CL.Z(ix) = IDX;
%     imageschs(DEM,CL)
%     
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017



X   = zeros(numel(varargin{1}.Z),nargin);
txt = cell(1,nargin);
txtname = cell(1,nargin);
for r=1:numel(varargin)
    
    if r>1
        validatealignment(varargin{r},varargin{1});
    end
    
    X(:,r) = double(varargin{r}.Z(:));
    txt{r} = inputname(r);
    txtname{r} = varargin{r}.name;
end

I  = ~any(isnan(X),2);
ix = (1:numel(varargin{1}.Z))';

X  = X(I,:);
ix = ix(I);
