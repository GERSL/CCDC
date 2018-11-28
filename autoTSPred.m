function outfity=autoTSPred(outfitx,fit_cft)
% Auto Trends and Seasonal Predict
% INPUTS:
% outfitx - Julian day [1; 2; 3];
% fit_cft - fitted coefficients;
% OUTPUTS:
% outfity - predicted reflectances [0.1; 0.2; 0.3];
% General model TSModel:
% f(x) =  a0 + b0*x + a1*cos(x*w) + b1*sin(x*w) 

% num_yrs = 365.25; % number of days per year
% w=2*pi/num_yrs; % anual cycle 
w = 2*pi/365.25;

outfity=[ones(size(outfitx)),outfitx,...% overall ref + trending
        cos(w*outfitx),sin(w*outfitx),...% add seasonality
        cos(2*w*outfitx),sin(2*w*outfitx),...% add bimodal seasonality
        cos(3*w*outfitx),sin(3*w*outfitx)]*fit_cft; % add trimodal seasonality
end
