function fit_cft=autoRobustFit(x,y,Nyr)
% INPUTS:
% x - Julian day [1; 2; 3];
% y - predicted reflectances [0.1; 0.2; 0.3];
% w - cofficient related with T=365d;
% yr - number of years ceil([end(day)-start(day)]/365);
% outfitx - fitted x values;

% OUTPUTS:
% outfity - fitted y values;
% fit_cft - fitted coefficients;
% General fitting model: f(x) =  a0 + a1*cos(x*w/N) + b1*sin(x*w/N) + ...
% a2*cos(x*w/1) + b2*sin(x*w/1)

% check num of clr obs before using
num_x = length(x); % number of clear pixels

% build X
X = zeros(num_x,4);
% annual cycle
w = 2*pi/365.25; % annual cycle
X(:,1) = cos(w*x);
X(:,2) = sin(w*x);
% Nyr year cycle
w = w/Nyr;
X(:,3) = cos(w*x);
X(:,4) = sin(w*x);

% Robust fitting
fit_cft = robustfit_cor(X,y);
end
