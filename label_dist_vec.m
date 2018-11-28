function [dist_type,dist_year,dist_doy]=label_dist_vec(curr_time,t_min,vec)
% This function is used to provde distubance year and disturbance type for
% each year based on NBR changes
% Version 1.0: (03/29/2018)
%% get disurbance pixel
% vec = obs - pred
% only provide disturbance map
% obs - pred
nir = vec(4);
vis = min(vec(1:3));
swir = min(vec([5,7]));

if nir > t_min && vis < -t_min && swir < -t_min 
    dist_type = 2;
else
    dist_type = 3; % land disturbance
end

dist_year = datevecmx(curr_time);
dist_year = dist_year(1);
dist_doy = curr_time - datenummx(dist_year,1,0);