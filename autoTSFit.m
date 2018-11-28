function [fit_cft,rmse,v_dif]=autoTSFit(x,y,df)
% Revisions: 
% v1.0 Using lasso for timeseries modeling (01/27/2013)
% Auto Trends and Seasonal Fit between breaks
% INPUTS:
% x - Julian day [1; 2; 3];
% y - predicted reflectances [0.1; 0.2; 0.3];
% df - degree of freedom (num_c)
%
% OUTPUTS:
% fit_cft - fitted coefficients;
% General model TSModel:
% f1(x) = a0 + b0*x (df = 2)
% f2(x) = f1(x) + a1*cos(x*w) + b1*sin(x*w) (df = 4)
% f3(x) = f2(x) + a2*cos(x*2w) + b2*sin(x*2w) (df = 6)
% f4(x) = f3(x) + a3*cos(x*3w) + b3*sin(x*3w) (df = 8)

n=length(x); % number of clear pixels
% num_yrs = 365.25; % number of days per year
w=2*pi/365.25; % num_yrs; % anual cycle
% fit coefs
fit_cft = zeros(8,1);

% global lamda_lasso

%% LASSO Fit
% build X
X = zeros(n,df-1);
X(:,1) = x;

if df >= 4
    X(:,2)=cos(w*x);
    X(:,3)=sin(w*x);
end

if df >= 6
    X(:,4)=cos(2*w*x);
    X(:,5)=sin(2*w*x);
end

if df >= 8
    X(:,6)=cos(3*w*x);
    X(:,7)=sin(3*w*x);
end

% lasso fit with lambda = 20
fit = glmnet_fast(X,y,glmnetSetL(20));  

% curr_cft=[fit.a0;fit.beta];
fit_cft(1:df) = [fit.a0;fit.beta]; % curr_cft;

% %% OLS Fit
% 
% % build X
% X = zeros(n,df);
% X(:,1) = ones(n,1);
% X(:,2) = x;
% 
% if df >= 4
%     X(:,3)=cos(w*x);
%     X(:,4)=sin(w*x);
% end
% 
% if df >= 6
%     X(:,5)=cos(2*w*x);
%     X(:,6)=sin(2*w*x);
% end
% 
% if df >= 8
%     X(:,7)=cos(3*w*x);
%     X(:,8)=sin(3*w*x);
% end 
% 
% % curr_cft=[fit.a0;fit.beta];
% fit_cft(1:df) = X\y; % curr_cft;

yhat=autoTSPred(x,fit_cft);
% rmse=median(abs(y-yhat));
v_dif = y-yhat;
rmse=norm(v_dif)/sqrt(n-df);
% f(x) = a0 + b0*x + a1*cos(x*w) + b1*sin(x*w) (df = 4)
end
