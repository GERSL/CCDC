function [dist_type,dist_year,dist_doy]=label_dist_type(curr_cft,curr_time,t_c,vec,next_cft)
% This function is used to provde distubance year and disturbance type 
% 1 => regrowth break
% 2 => aforestation break
% 3 => land disturbance
% 
% Version 1.0: (Zhe Zhu 10/30/2018)
%% get disurbance pixel
% vec = obs - pred
% only provide disturbance map
% obs - pred
nir = vec(4);
c_nir = curr_cft(2,4);
n_nir = next_cft(2,4);

vis = vec(3);
c_vis = curr_cft(2,3);
n_vis = next_cft(2,3);

swir = vec(5);
c_swir = curr_cft(2,5);
n_swir = next_cft(2,5);

if nir > t_c && vis < -t_c && swir < -t_c
    if c_nir > abs(n_nir) && c_vis < -abs(n_vis) && c_swir < -abs(n_swir)
        dist_type = 2; % aforestation
    else
        dist_type = 1; % regrowth
    end
else
    dist_type = 3; % land disturbance
end

dist_year = datevecmx(curr_time);
dist_year = dist_year(1);
dist_doy = curr_time - datenummx(dist_year,1,0);