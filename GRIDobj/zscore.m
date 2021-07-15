function X = zscore(X)

%ZSCORE standardized z-scores for GRIDobj
%
% Syntax
%
%     Z = zscore(X)
%
% Description
%
%     zscore returns the z-score for each element of GRIDobj X such that 
%     all values of X are centered to have mean 0 and scaled to have 
%     standard deviation 1.
%
% Input arguments
%
%     X     GRIDobj
%
% Output arguments
%
%     Y     GRIDobj
% 
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     Z = zscore(DEM);
%     imagesc(Z); colorbar
% 
% See also: GRIDobj, zscore
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 15. March, 2016



I = ~isnan(X.Z);
X.Z(I) = (X.Z(I(:)) - mean(X.Z(I(:))))./std(X.Z(I(:)));