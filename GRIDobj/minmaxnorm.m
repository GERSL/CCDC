function A = minmaxnorm(A,prcclip)

%MINMAXNORM min-max normalization with optional percent clipping
%
% Syntax
%
%     A = minmaxnorm(A)
%     A = minmaxnorm(A,prc)
%
% Description
%
%     minmaxnorm normalizes GRIDobj A to a range between 0 and 1.
%     Optionally, extremes of the data can be removed by percent clipping
%     such that values higher (lower) the 1-prc (prc) percentile are set 
%     to 1 (0).
%
% Input arguments
%
%     A     GRIDobj
%     prc   percentile, scalar value between 0 and 100.
%

if nargin == 1
    minz = min(A);
    maxz = max(A);
    A = (A - minz)/(maxz-minz);
else
    Z = A.Z;
    
    qclip = prcclip/100;
    [n,edges] = histcounts(Z(~isnan(Z(:))),'Normalization','cdf');
    lval = edges(find(n>=qclip,1,'first'));
    uval = edges(find(n<(1-qclip),1,'last'));
    if lval == uval
        warning('TopoToolbox:minmaxnorm','percent clip returns flat matrix');
        Z(:,:) = lval;
    else
        Z = max(Z,lval);
        Z = min(Z,uval);
    end
    
    Z = (Z-lval)/(uval-lval);
    A.Z = Z;
    
end

