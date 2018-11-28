function b = statrobustfit_cor(X,y,wfun,tune)
% STATROBUSTFIT Calculation function for ROBUSTFIT
% Corrected 1.0 version (Zhe 04/06/2013)

[n,p] = size(X);
X = [ones(n,1) X];
p = p+1;

% Find the least squares solution.
[Q,R,perm] = qr(X,0);
tol = abs(R(1)) * max(n,p) * eps(class(R));
xrank = sum(abs(diag(R)) > tol);
if xrank==p
    b(perm,:) = R \ (Q'*y);
else
    % Use only the non-degenerate parts of R and Q, but don't reduce
    % R because it is returned in stats and is expected to be of
    % full size.
    b(perm,:) = [R(1:xrank,1:xrank) \ (Q(:,1:xrank)'*y); zeros(p-xrank,1)];
    perm = perm(1:xrank);
end
b0 = zeros(size(b));

% Adjust residuals using leverage, as advised by DuMouchel & O'Brien
E = X(:,perm)/R(1:xrank,1:xrank);
h = min(.9999, sum(E.*E,2));
adjfactor = 1 ./ sqrt(1-h);

% If we get a perfect or near perfect fit, the whole idea of finding
% outliers by comparing them to the residual standard deviation becomes
% difficult.  We'll deal with that by never allowing our estimate of the
% standard deviation of the error term to get below a value that is a small
% fraction of the standard deviation of the raw response values.
tiny_s = 1e-6 * std(y);
if tiny_s == 0
    tiny_s = 1;
end

% Perform iteratively reweighted least squares to get coefficient estimates
D = sqrt(eps(class(X)));
iter = 1;
iterlim = 5;% default value is 50 (less iteration improves efficiency)
wxrank = xrank;    % rank of weighted version of x
while((iter==0) || any(abs(b-b0) > D*max(abs(b),abs(b0))))
    iter = iter+1;
    if (iter>iterlim)
        break;
    end
    
    % Compute residuals from previous fit, then compute scale estimate
    r = y - X*b;
    radj = r .* adjfactor;
    mad_s = madsigma(radj,wxrank);% oringinally "s"
    
    % Compute new weights from these residuals, then re-fit
    w = feval(wfun, radj/(max(mad_s,tiny_s)*tune));% oringinally "s"
    b0 = b;
    [b(perm),wxrank] = wfit(y,X(:,perm),w);
end

% -----------------------------
function [b,r] = wfit(y,x,w)
%WFIT    weighted least squares fit

% Create weighted x and y
n = size(x,2);
sw = sqrt(w);
yw = y .* sw;
xw = x .* sw(:,ones(1,n));

% Computed weighted least squares results
[b,r] = linsolve(xw,yw,struct('RECT',true));

% -----------------------------
function s = madsigma(r,p)
%MADSIGMA    Compute sigma estimate using MAD of residuals from 0
rs = sort(abs(r));
s = median(rs(max(1,p):end)) / 0.6745;

