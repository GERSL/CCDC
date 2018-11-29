function rec_cg = TrendSeasonalFit_v12_30Line(dir_l,n_rst,ncols,nrows,T_cg,Tmax_cg,conse,num_c,nbands,B_detect)
% CCDC 12.30 version - Zhe Zhu, EROS USGS
% It is based on 7 bands fitting for Iterative Seasonal, Linear, and Break Models
% This function works for analyzing one line of time series pixel
%
%% Revisions: $ Date: 11/20/2015 $ Copyright: Zhe Zhu
%  Version 12.30  Fixed a bug for pixels without minimum observations (11/20/2015)
%  Version 12.29  Modified fit for perennial snow and Fmask failed pixels (09/20/2015)
%  Version 12.29  Do not fit disturbed time period (09/18/2015)
%  Version 12.28  Fixed a bug for missing values in land cover maps (09/16/2015)
%  Version 12.27  Fixed bugs for persistent snow and falied Fmask pixels (06/17/2015)
%  Version 12.26  Connected time for all models (05/28/2015)
%  Version 12.25  Bug fixed in snow percent (05/19/2015)
%  Version 12.24  Change T_const in Tmask (03/31/2015)
%  Version 12.23  Update iteratively before 24 observations (03/22/2015)
%  Version 12.22  Adjust mini RMSE based on temporal variability (03/22/2015)
%  Version 12.21  Add more categories and update i_start in the end (03/14/2015)
%  Version 12.20  Convert BT to from K to C before analysis (03/12/2015)
%  Version 12.19  Fit for permanent snow if is more than 75% (03/12/2015)
%  Version 12.18  No change detection if clear observation less than 25% (03/12/2015)
%  Version 12.17  Use median value for very simple model & change magnitude (02/24/2015)
%  Version 12.16  Finding changes in all water pixels (02/24/2015)
%  Version 12.15  Use the original multitemporal cloud mask (02/15/2015)
%  Version 12.14  Do not need example_img in images folder (02/09/2015)
%  Version 12.13: More infromation in "category" (11/10/2014)
%  This version (12.13) is used for the third round of the LCMS project.
%  Command: TrendSeasonalFit_v12Plot(N_row,N_col,min=0.5,T_cg=0.99,n_times=3,conse=6,B_detect=2:6)
%  Version 12.12: Fit for pixels where Fmask fails (11/09/2014)
%  Version 12.11: Bug fixed in num_fc (11/09/2014)
%  Version 12.10: Better multietmporal cloud detection at the beginning (11/06/2014)
%  Version 12.9:  Detect change only for land pixel (water/snow speical case) (10/31/2014)
%  Version 12.8:  Speed up by reducing time for RMSE and model computing (10/17/2014)
%  Version 12.7:  mini rmse should be larger than 10% of the mean (10/13/2014)
%  Version 12.6:  Fit model again when there are a 33.3% more data (10/08/2014)
%  Version 12.5:  Use subset of bands (2-6) for detecting surface change (10/01/2014)
%  Version 12.4:  Only apply multitemporal cloud masking during model initialization (09/29/2014)
%  Version 12.3:  Use subset of bands (3-5) to balance change in diferent dimensions (09/01/2014)
%  This version (12.3) is used for the second round of the LCMS project.
%  Command: TrendSeasonalFit_v12Plot(N_row,N_col,min=1,T_cg=0.99,n_times=3,conse=5,B_detect=3:6)
%  Version 12.2:  Bug fixed in model intialization (08/14/2014)
%  Version 12.1:  Use subset of bands (3-6) to avoid atmosphere influences (08/04/2014)
%% Version 12.0   Detecting change based on probability (07/19/2014)
%  Version 11.6:  No need to change folder name & faster in speed (by Christ Holden 06/06/2014)
%  Version 11.5:  Improved calculation of temporally adjusted RMSE (04/23/2014)
%  Version 11.4:  Revise "rec_cg.category" to better seperate different fit processes (04/01/2014)
%  This version (11.4) is used for generating synthetic data for ACRE project and
%  detecting change for LCMS project.
%  Command: TrendSeasonalFit_v11Plot(N_row,N_col,min=1,T_cg=2,n_times=3,conse=6,B_detect=1:6)
%  Version 11.3:  Add "rec_cg.magnitude" as indicator of change magnitude (04/01/2014)
%  Version 11.2:  Change very simple fit with mean value for start and end of timeseries (04/01/2014)
%  Version 11.1:  Do not need metadata in the image folder to run CCDC (03/25/2014)
%% Version 11.0:  Use change vector magnitude as threshold for detecting change (03/25/2014)
%  Version 10.13: Use subset of bands (1-6) to avoid atmosphere influences (01/31/2014)
%  Version 10.12: More accurate number of days per year "num_yrs" (01/30/2014)
%  Version 10.11: RMSE updates with time series fit (01/26/2014)
%  Version 10.10: Update temperature extreme in recent studies (01/16/2014)
%  Version 10.9:  Find break in max value in any of the band (01/08/2014)
%  Version 10.8:  Add very simple fit with median value for start and end of timeseries (10/21/2013)
%  This version (10.8) is used for generating synthetic data for the LCMS project.
%  Command: TrendSeasonalFit_v10Plot('stack',N_row,N_col,mini=0.5,T_cg=3,n_times=3,conse=6,B_detect=2:6)
%  Version 10.7:  Better multitemporal cloud detection (10/19/2013)
%  Version 10.6:  Add "Tmax_cg" for last step noise removal (10/18/2013)
%  Version 10.5:  Use subset of bands (2-6) to avoid atmosphere influences (10/18/2013)
%  Version 10.4:  Let dynamic fitting for pixels at the beginning (09/23/2013)
%  Version 10.3:  Able to detect change at the verying beginning (09/06/2013)
%  Version 10.2:  Add mini years "mini_yrs" in model intialization (09/03/2013)
%  Version 10.1:  Reduce time for calcuating "v_dif" (09/02/2013)
%% Version 10.0:  Fit for beginning and end of the time series (08/31/2013)
%  Version 9.9:   Only fit more than 50% of Landat images overlap area (08/28/2013)
%  Version 9.8:   Force model fit for persistent snow pixels (08/27/2013)
%  Version 9.7:   Add "rec_cg.category" as indicator of fitting procudure (08/20/2013)
%                 Add rec_cg.change_prob as indicator of change probability (08/20/2013)
%                 Add rec_cg.num_obs ad indicator of number of observations (08/20/2013)
%  Version 9.6:   Remove mininum rmse "mini" and minimum years "mini_yrs" (08/16/2013)
%  Version 9.5:   Model gets more coefficients with more observations (08/16/2013)
%  Version 9.4:   Bug fixed in calculating temporally adjusted rmse (08/01/2013)
%  Version 9.3:   Fit curve again after one year (03/28/2013)
%  This version (9.3) is used for mapping land cover for the IDS project.
%  Command: TrendSeasonalFit_v9Plot('stack',N_row,N_col,T_cg=2,n_times=3,conse=4)
%  Version 9.2:   Use "mini = T_const/T_cg" for small rmse cases (03/26/2013)
%  Version 9.1:   Remove out of range pixels before time series analysis (02/09/2013)
%% Version 9.0:   Using 8 coefficients and lasso fit (02/01/2013)
%  Version 8.4:   Use "max v_slope" instead of "average v_slope" (01/16/2013)
%  Version 8.3:   Start initialization when "time_span" > 1 year (01/16/2013)
%  Version 8.2:   Bug fixed in not fitting models at the begining (01/16/2013)
%  Version 8.1:   Bug fixed in counting "i" and "i_span"(01/13/2013)
%% Version 8.0:   Temporally changing RMSE (01/09/2013)
%% Version 7.3:   Continuous Change Detection and Classification (CCDC) (07/11/2012)
%  This version (7.3) is explained by Zhu, Z. & Woodcock, C.E., Continuous Change
%  Detection and Classification (CCDC) of land cover using all available
%  Landsat data, Remote Sensing of Environment (2014).
%  Command: TrendSeasonalFit_v7Plot('stack',N_row,N_col,T_cg=3,n_times=3,conse=3)
%% Version 1.0:   Continous Monitoring of Forest Disturbance Algorithm (CMFDA) (07/13/2010)
%  This version (1.0) is explained by Zhu, Z., Woodcock, C.E., Olofsson, P.,
%  Continuous monitoring of forest disturbance using all available Landsat
%  data, Remote Sensing of Environment (2012).
%
%% Inputs:
% stk_n='stack'; stack image name
% ncols = 8021; % number of pixels processed per line
% nrows=1; % the nrowsth lines
% for example    1 2 3 4 5
%                6 7 8 9 10
%
%% Outputs:
%
% rec_cg RECord information about all curves between ChanGes
% rec_cg(i).t_start record the start of the ith curve fitting (julian_date)
% rec_cg(i).t_end record the end of the ith curve fitting (julian_date)
% rec_cg(i).t_break record the first observed break time (julian_date)
% rec_cg(i).coefs record the coefficients of the ith curve
% rec_cg(i).pos record the position of the ith pixel (pixel id)
% rec_cg(i).magnitude record the change vector of all spectral bands
% rec_cg(i).category record what fitting procudure and model is used
% cateogry category 5x: persistent snow    4x: Fmask fails
% cateogry category 3x: modified fit       2x: end fit
% category category 1x: start fit           x: normal procedure
% cateogry category x1: mean value         x4: simple model
% category category x6: advanced model     x8: full model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  defining variables
%% Constants
% maximum number of coefficient required
% 2 for tri-modal; 2 for bi-modal; 2 for seasonality; 2 for linear;
min_num_c = 4;
mid_num_c = 6;
max_num_c = 8;
% number of clear observation / number of coefficients
n_times = 3;
% initialize NUM of Functional Curves
num_fc = 0;
% number of days per year
num_yrs = 365.25;
% number of bytes: int16
num_byte = 2;
% Band for multitemporal cloud/snow detection (Green)
num_B1 = 2;
% Band for multitemporal shadow/snow shadow detection (SWIR)
num_B2 = 5;
% Threshold for cloud, shadow, and snow detection.
T_const = 3.89;
% minimum year for model intialization
mini_yrs = 1;
% no change detection for permanent snow pixels
t_sn = 0.75;
% Fmask fails threshold
t_clr = 0.25;

% get num of total folders start with "L"
imf = dir(fullfile(dir_l,'L*')); % folder names
% filter for Landsat folders
imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
imf = [imf{:}];
imf = vertcat(imf{:});
% sort according to yeardoy
yeardoy = str2num(imf(:, 10:16));
[~, sort_order] = sort(yeardoy);
imf = imf(sort_order, :);
% number of folders start with "L"
num_t = size(imf,1);

% initialize the struct data of RECording of ChanGe (rec_cg)
rec_cg = struct('t_start',[],'t_end',[],'t_break',[],'coefs',[],'rmse',[],...
    'pos',[],'change_prob',[],'num_obs',[],'category',[],'magnitude',[]);

% % mask for study area (1 fit, 0 no fit)
% fit_mask = enviread('GZ_Mask');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Get ready for Xs & Ys
%% Read in Xs & Ysq
% transforming to serial date number (0000 year)
sdate = zeros(num_t,1); % Xs
line_t = zeros(num_t,nbands*ncols); %Ys

for i=1:num_t
    im_dir = dir(fullfile(dir_l,imf(i, :)));
    im = '';
    for f = 1:size(im_dir, 1)
        % use regular expression to match:
        %   'L(\w*)'    Any word begining with L that has any following chars
        %   stk_n       includes stack name somewhere after L
        %   '$'           ends with the stack name (e.g., no .hdr, .aux.xml)
        if regexp(im_dir(f).name, ['L(\w*)', 'stack', '$']) == 1
            im = fullfile(dir_l,imf(i, :), im_dir(f).name);
            break
        end
    end
    % Check to make sure we found something
    if strcmp(im, '')
        error('Could not find stack image for directory %s\n', imf(i));
    end
    % Find date for folder imf(i)
    yr = str2num(imf(i, 10:13));
    doy = str2num(imf(i, 14:16));
    sdate(i) = datenum(yr, 1, 0) + doy;
    dummy_name = im;
    fid_t = fopen(dummy_name,'r'); % get file ids
    fseek(fid_t,num_byte*(nrows-1)*ncols*nbands,'bof');
    line_t(i,:) = fread(fid_t,nbands*ncols,'int16=>double','ieee-le'); % get Ys
end
fclose('all'); % close all files

for i_ids = 1:ncols
    % mask data
    line_m = line_t(:,nbands*i_ids);
    
    % Only run CCDC for places where more than 50% of images has data
    idexist = line_m < 255;
    overlap_pct = sum(idexist)/num_t;
    if overlap_pct < 0.5
        continue;
    end
    
    % convert Kelvin to Celsius
    line_t(:,nbands*(i_ids-1)+7) = line_t(:,nbands*(i_ids-1)+7)*10 - 27315;
    
    % clear pixel should have reflectance between 0 and 1
    % brightness temperature should between -93.2 to 70.7 celsius degree
    idrange = line_t(:,nbands*(i_ids-1)+1)>0&line_t(:,nbands*(i_ids-1)+1)<10000&...
        line_t(:,nbands*(i_ids-1)+2)>0&line_t(:,nbands*(i_ids-1)+2)<10000&...
        line_t(:,nbands*(i_ids-1)+3)>0&line_t(:,nbands*(i_ids-1)+3)<10000&...
        line_t(:,nbands*(i_ids-1)+4)>0&line_t(:,nbands*(i_ids-1)+4)<10000&...
        line_t(:,nbands*(i_ids-1)+5)>0&line_t(:,nbands*(i_ids-1)+5)<10000&...
        line_t(:,nbands*(i_ids-1)+6)>0&line_t(:,nbands*(i_ids-1)+6)<10000&...
        line_t(:,nbands*(i_ids-1)+7)>-9320&line_t(:,nbands*(i_ids-1)+7)<7070;
    
    % # of clear observatons
    idclr = line_m < 2;
    % # of all available observations
    idall = line_m < 255;
    % clear observation percentage
    clr_pct = sum(idclr)/sum(idall);
    % snow pixels
    idsn = line_m == 3;
    % percent of snow observations
    sn_pct = sum(idsn)/(sum(idclr)+sum(idsn)+0.01);
    
    % not enough clear observations for change detection
    if clr_pct < t_clr
        % permanent snow pixels
        if sn_pct > t_sn
            % snow observations are "good" now
            idgood = idsn|idclr;
            % number of snow pixel within range
            n_sn = sum(idgood);
            
            if n_sn < n_times*min_num_c % not enough snow pixels
                continue
            else
                % Xs & Ys for computation
                clrx = sdate(idgood);
                % bands 1-5,7,6
                clry = line_t(idgood,(nbands*(i_ids-1)+1):(nbands*(i_ids-1)+nbands-1));
                
                % the first observation for TSFit
                i_start = 1;
                % identified and move on for the next curve
                num_fc = num_fc + 1; % NUM of Fitted Curves (num_fc)
                
                % defining computed variables
                fit_cft = zeros(max_num_c,nbands-1);
                % rmse for each band
                rmse = zeros(nbands-1,1);
                % snow qa = 50
                qa = 50;
                
                for i_B=1:nbands-1
                    if i_B ~= nbands-1 % treat saturated and unsaturated pixels differently
                        idgood = clry(:,i_B) < 10000; % saturate if ref > 1 or NBR NDVR > 1
                        i_span = sum(idgood);
                        if i_span < min_num_c*n_times % fill value for frequently saturated snow pixels
                            fit_cft(1,i_B) = 10000; % give constant value
                        else % fit for enough unsaturat snow pixels
                            [fit_cft(:,i_B),rmse(i_B)]=autoTSFit(clrx(idgood),clry(idgood,i_B),min_num_c);
                        end
                    else % fit for temperature band
                        idgood = clry(:,i_B)>-9320&clry(:,i_B)<7070;
                        [fit_cft(:,i_B),rmse(i_B)]=autoTSFit(clrx(idgood),clry(idgood,i_B),min_num_c);
                    end
                end
                
                % updating information at each iteration
                % record time of curve start
                rec_cg(num_fc).t_start=clrx(i_start);
                % record time of curve end
                rec_cg(num_fc).t_end=clrx(end);
                % record break time
                rec_cg(num_fc).t_break = 0; % no break at the moment
                % record postion of the pixel
                rec_cg(num_fc).pos = (nrows-1)*ncols+i_ids;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = n_sn;
                % record fit category
                rec_cg(num_fc).category = qa + min_num_c;
                % record change magnitude
                rec_cg(num_fc).magnitude = zeros(1,nbands-1);
            end
        else % no change detection for clear observations
            % within physical range pixels
            idgood = idrange;
            
            % Xs & Ys for computation
            clrx = sdate(idgood);
            % bands 1-5,7,6
            clry = line_t(idgood,(nbands*(i_ids-1)+1):(nbands*(i_ids-1)+nbands-1));
            
            idclr = clry(:,num_B1) < median(clry(:,num_B1)) + 400;
            n_clr = sum(idclr);
            
            if n_clr < n_times*min_num_c % not enough clear pixels
                continue
            else
                % Xs & Ys for computation
                clrx = clrx(idclr);
                clry = clry(idclr,:);
                
                % the first observation for TSFit
                i_start = 1;
                % identified and move on for the next curve
                num_fc = num_fc + 1; % NUM of Fitted Curves (num_fc)
                
                % defining computed variables
                fit_cft = zeros(max_num_c,nbands-1);
                % rmse for each band
                rmse = zeros(nbands-1,1);
                % Fmask fail qa = 40
                qa = 40;
                
                for i_B = 1:nbands-1
                    % fit basic model for all within range snow pixels
                    [fit_cft(:,i_B),rmse(i_B)] = autoTSFit(clrx,clry(:,i_B),min_num_c);
                end
                
                % record time of curve start
                rec_cg(num_fc).t_start = clrx(i_start);
                % record time of curve end
                rec_cg(num_fc).t_end = clrx(end);
                % record break time
                rec_cg(num_fc).t_break = 0;
                % record postion of the pixel
                rec_cg(num_fc).pos = (nrows-1)*ncols+i_ids;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = length(clrx);
                % record fit category
                rec_cg(num_fc).category = qa + min_num_c;
                % record change magnitude
                rec_cg(num_fc).magnitude = zeros(1,nbands-1);
            end
        end
    else % normal CCDC procedure
        
        % clear and within physical range pixels
        idgood = idclr & idrange;
        
        % Xs & Ys for computation
        clrx = sdate(idgood);
        % bands 1-5,7,6
        clry = line_t(idgood,(nbands*(i_ids-1)+1):(nbands*(i_ids-1)+nbands-1));
        
        % caculate median variogram
        var_clry = clry(2:end,:)-clry(1:end-1,:);
        adj_rmse = median(abs(var_clry),1);
        
        % start with the miminum requirement of clear obs
        i = n_times*min_num_c;
        
        % initializing variables
        % the first observation for TSFit
        i_start = 1;
        % record the start of the model initialization (0=>initial;1=>done)
        BL_train = 0;
        % identified and move on for the next curve
        num_fc = num_fc+1; % NUM of Fitted Curves (num_fc)
        % record the num_fc at the beginning of each pixel
        rec_fc = num_fc;
        
        % while loop - process till the last clear observation - conse
        while i<= length(clrx)-conse
            % span of "i"
            i_span = i-i_start+1;
            % span of time (num of years)
            time_span=(clrx(i)-clrx(i_start))/num_yrs;
            
            % basic requrirements: 1) enough observations; 2) enough time
            if i_span >= n_times*min_num_c && time_span >= mini_yrs
                % initializing model
                if BL_train == 0
                    % Tmask: noise removal (good => 0 & noise => 1)
                    blIDs = autoTmask(clrx(i_start:i+conse),clry(i_start:i+conse,[num_B1,num_B2]),...
                        (clrx(i+conse)-clrx(i_start))/num_yrs,adj_rmse(num_B1),adj_rmse(num_B2),T_const);
                    
                    % IDs to be removed
                    IDs = i_start:i+conse;
                    rmIDs = IDs(blIDs(1:end-conse) == 1);
                    
                    % update i_span after noise removal
                    i_span = sum(~blIDs(1:end-conse));
                    
                    % check if there is enough observation
                    if i_span < n_times*min_num_c
                        % move forward to the i+1th clear observation
                        i = i+1;
                        % not enough clear observations
                        continue;
                        % check if there is enough time
                    else
                        % copy x & y
                        cpx = clrx;
                        cpy = clry;
                        
                        % remove noise pixels between i_start & i
                        cpx(rmIDs) = [];
                        cpy(rmIDs,:) = [];
                        
                        % record i before noise removal
                        % This is very important as if model is not initialized
                        % the multitemporal masking shall be done again instead
                        % of removing outliers in every masking
                        i_rec = i;
                        
                        % update i afer noise removal (i_start stays the same)
                        i = i_start+i_span-1;
                        % update span of time (num of years)
                        time_span = (cpx(i)-cpx(i_start))/num_yrs;
                        
                        % check if there is enough time
                        if time_span < mini_yrs
                            % keep the original i
                            i = i_rec;
                            % move forward to the i+1th clear observation
                            i = i+1;
                            % not enough time
                            continue;
                            % Step 2: model fitting
                        else
                            % remove noise
                            clrx = cpx;
                            clry = cpy;
                            
                            % Step 2: model fitting
                            % initialize model testing variables
                            % defining computed variables
                            fit_cft = zeros(max_num_c,nbands-1);
                            % rmse for each band
                            rmse = zeros(nbands-1,1);
                            % value of differnce
                            v_dif = zeros(nbands-1,1);
                            % record the diference in all bands
                            rec_v_dif = zeros(i-i_start+1,nbands-1);
                            
                            for i_B = 1:nbands-1
                                % initial model fit
                                [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                    autoTSFit(clrx(i_start:i),clry(i_start:i,i_B),min_num_c);
                            end
                            
                            % normalized to z-score
                            for i_B = B_detect
                                % minimum rmse
                                mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                                
                                % compare the first clear obs
                                v_start = rec_v_dif(1,i_B)/mini_rmse;
                                % compare the last clear observation
                                v_end = rec_v_dif(end,i_B)/mini_rmse;
                                % anormalized slope values
                                v_slope = fit_cft(2,i_B)*(clrx(i)-clrx(i_start))/mini_rmse;
                                
                                % differece in model intialization
                                v_dif(i_B) = abs(v_slope) + abs(v_start) + abs(v_end);
                            end
                            v_dif = norm(v_dif(B_detect))^2;
                            
                            % find stable start for each curve
                            if v_dif > T_cg
                                % start from next clear obs
                                i_start = i_start + 1;
                                % move forward to the i+1th clear observation
                                i = i + 1;
                                % keep all data and move to the next obs
                                continue;
                            else
                                % model ready!
                                BL_train = 1;
                                % count difference of i for each iteration
                                i_count = 0;
                                
                                % find the previous break point
                                if num_fc == rec_fc
                                    % first curve
                                    i_break = 1;
                                else
                                    % after the first curve
                                    i_break = find(clrx >= rec_cg(num_fc-1).t_break);
                                    i_break = i_break(1);
                                end
                                
                                if i_start > i_break
                                    % model fit at the beginning of the time series
                                    for i_ini = i_start-1:-1:i_break
                                        if i_start - i_break < conse
                                            ini_conse = i_start - i_break;
                                        else
                                            ini_conse = conse;
                                        end
                                        % value of difference for conse obs
                                        v_dif = zeros(ini_conse,nbands-1);
                                        % record the magnitude of change
                                        v_dif_mag = v_dif;
                                        % chagne vector magnitude
                                        vec_mag = zeros(ini_conse,1);
                                        
                                        for i_conse = 1:ini_conse
                                            for i_B = 1:nbands-1
                                                % absolute difference
                                                v_dif_mag(i_conse,i_B) = clry(i_ini-i_conse+1,i_B)-autoTSPred(clrx(i_ini-i_conse+1),fit_cft(:,i_B));
                                                % normalized to z-scores
                                                if sum(i_B == B_detect)
                                                    % minimum rmse
                                                    mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                                                    
                                                    % z-scores
                                                    v_dif(i_conse,i_B) = v_dif_mag(i_conse,i_B)/mini_rmse;
                                                end
                                            end
                                            vec_mag(i_conse) = norm(v_dif(i_conse,B_detect))^2;
                                        end
                                        
                                        if min(vec_mag) > T_cg % change detected
                                            break
                                        elseif vec_mag(1) > Tmax_cg % false change
                                            % remove noise
                                            clrx(i_ini) = [];
                                            clry(i_ini,:) = [];
                                            i=i-1; % stay & check again after noise removal
                                        end
                                        
                                        % update new_i_start if i_ini is not a confirmed break
                                        i_start = i_ini;
                                    end
                                end
                                
                                if num_fc == rec_fc && i_start - i_break >= conse
                                    % defining computed variables
                                    fit_cft = zeros(max_num_c,nbands-1);
                                    % rmse for each band
                                    rmse = zeros(nbands-1,1);
                                    % start fit qa = 10
                                    qa = 10;
                                    
                                    for i_B=1:nbands-1
                                        [fit_cft(:,i_B),rmse(i_B)] = ...
                                            autoTSFit(clrx(i_break:i_start-1),clry(i_break:i_start-1,i_B),min_num_c);
                                    end
                                    
                                    % record time of curve end
                                    rec_cg(num_fc).t_end = clrx(i_start-1);
                                    % record postion of the pixel
                                    rec_cg(num_fc).pos = (nrows-1)*ncols + i_ids;
                                    % record fitted coefficients
                                    rec_cg(num_fc).coefs = fit_cft;
                                    % record rmse of the pixel
                                    rec_cg(num_fc).rmse = rmse;
                                    % record break time
                                    rec_cg(num_fc).t_break = clrx(i_start);
                                    % record change probability
                                    rec_cg(num_fc).change_prob = 1;
                                    % record time of curve start
                                    rec_cg(num_fc).t_start = clrx(1);
                                    % record fit category
                                    rec_cg(num_fc).category = qa + min_num_c;
                                    % record number of observations
                                    rec_cg(num_fc).num_obs = i_start - i_break;
                                    % record change magnitude
                                    rec_cg(num_fc).magnitude = - median(v_dif_mag,1);
                                    
                                    % identified and move on for the next functional curve
                                    num_fc = num_fc + 1;
                                end
                            end
                        end
                    end
                end
                
                % continuous monitoring started!!!
                if BL_train == 1
                    % all IDs
                    IDs = i_start:i;
                    % span of "i"
                    i_span = i-i_start+1;
                    
                    % determine the time series model
                    update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c);
                    
                    % initial model fit when there are not many obs
                    if  i_count == 0 || i_span <= max_num_c*n_times
                        % update i_count at each interation
                        i_count = clrx(i)-clrx(i_start);
                        
                        % defining computed variables
                        fit_cft = zeros(max_num_c,nbands-1);
                        % rmse for each band
                        rmse = zeros(nbands-1,1);
                        % record the diference in all bands
                        rec_v_dif = zeros(length(IDs),nbands-1);
                        % normal fit qa = 0
                        qa = 0;
                        
                        for i_B=1:nbands-1
                            [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                autoTSFit(clrx(IDs),clry(IDs,i_B),update_num_c);
                        end
                        
                        % updating information for the first iteration
                        % record time of curve start
                        rec_cg(num_fc).t_start = clrx(i_start);
                        % record time of curve end
                        rec_cg(num_fc).t_end = clrx(i);
                        % record break time
                        rec_cg(num_fc).t_break = 0; % no break at the moment
                        % record postion of the pixel
                        rec_cg(num_fc).pos = (nrows-1)*ncols+i_ids;
                        % record fitted coefficients
                        rec_cg(num_fc).coefs = fit_cft;
                        % record rmse of the pixel
                        rec_cg(num_fc).rmse = rmse;
                        % record change probability
                        rec_cg(num_fc).change_prob = 0;
                        % record number of observations
                        rec_cg(num_fc).num_obs = i-i_start+1;
                        % record fit category
                        rec_cg(num_fc).category = qa + update_num_c;
                        % record change magnitude
                        rec_cg(num_fc).magnitude = zeros(1,nbands-1);
                        
                        % detect change
                        % value of difference for conse obs
                        v_dif = zeros(conse,nbands-1);
                        % record the magnitude of change
                        v_dif_mag = v_dif;
                        vec_mag = zeros(conse,1);
                        
                        for i_conse = 1:conse
                            for i_B = 1:nbands-1
                                % absolute difference
                                v_dif_mag(i_conse,i_B) = clry(i+i_conse,i_B)-autoTSPred(clrx(i+i_conse),fit_cft(:,i_B));
                                % normalized to z-scores
                                if sum(i_B == B_detect)
                                    % minimum rmse
                                    mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                                    
                                    % z-scores
                                    v_dif(i_conse,i_B) = v_dif_mag(i_conse,i_B)/mini_rmse;
                                end
                            end
                            vec_mag(i_conse) = norm(v_dif(i_conse,B_detect))^2;
                        end
                        % IDs that haven't updated
                        IDsOld = IDs;
                    else
                        if clrx(i)-clrx(i_start) >= 1.33*i_count
                            % update i_count at each interation
                            i_count = clrx(i)-clrx(i_start);
                            
                            % defining computed variables
                            fit_cft = zeros(max_num_c,nbands-1);
                            % rmse for each band
                            rmse = zeros(nbands-1,1);
                            % record the diference in all bands
                            rec_v_dif = zeros(length(IDs),nbands-1);
                            % normal fit qa = 0
                            qa = 0;
                            
                            for i_B=1:nbands-1
                                [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                    autoTSFit(clrx(IDs),clry(IDs,i_B),update_num_c);
                            end
                            
                            % record fitted coefficients
                            rec_cg(num_fc).coefs = fit_cft;
                            % record rmse of the pixel
                            rec_cg(num_fc).rmse = rmse;
                            % record number of observations
                            rec_cg(num_fc).num_obs = i-i_start+1;
                            % record fit category
                            rec_cg(num_fc).category = qa + update_num_c;
                            
                            % IDs that haven't updated
                            IDsOld = IDs;
                        end
                        
                        % record time of curve end
                        rec_cg(num_fc).t_end = clrx(i);
                        
                        % use fixed number for RMSE computing
                        n_rmse = n_times*rec_cg(num_fc).category;
                        tmpcg = zeros(nbands-1,1);
                        % better days counting for RMSE calculating
                        % relative days distance
                        d_rt = clrx(IDsOld) - clrx(i+conse);
                        d_yr = abs(round(d_rt/num_yrs)*num_yrs-d_rt);
                        
                        [~,sorted_indx] = sort(d_yr);
                        sorted_indx = sorted_indx(1:n_rmse);
                        
                        for i_B = B_detect
                            % temporally changing RMSE
                            tmpcg_rmse(i_B) = norm(rec_v_dif(IDsOld(sorted_indx)-IDsOld(1)+1,i_B))/...
                                sqrt(n_rmse-rec_cg(num_fc).category);
                        end
                        
                        % move the ith col to i-1th col
                        v_dif(1:conse-1,:) = v_dif(2:conse,:);
                        % only compute the difference of last consecutive obs
                        v_dif(conse,:) = 0;
                        % move the ith col to i-1th col
                        v_dif_mag(1:conse-1,:) = v_dif_mag(2:conse,:);
                        % record the magnitude of change of the last conse obs
                        v_dif_mag(conse,:) = 0;
                        % move the ith col to i-1th col
                        vec_mag(1:conse-1) = vec_mag(2:conse);
                        % change vector magnitude
                        vec_mag(conse) = 0;
                        
                        for i_B = 1:nbands-1
                            % absolute difference
                            v_dif_mag(conse,i_B) = clry(i+conse,i_B)-autoTSPred(clrx(i+conse),fit_cft(:,i_B));
                            % normalized to z-scores
                            if sum(i_B == B_detect)
                                % minimum rmse
                                mini_rmse = max(adj_rmse(i_B),tmpcg_rmse(i_B));
                                
                                % z-scores
                                v_dif(conse,i_B) = v_dif_mag(conse,i_B)/mini_rmse;
                            end
                        end
                        vec_mag(conse) = norm(v_dif(end,B_detect))^2;
                    end
                    
                    % change detection
                    if min(vec_mag) > T_cg % change detected
                        % record break time
                        rec_cg(num_fc).t_break = clrx(i+1);
                        % record change probability
                        rec_cg(num_fc).change_prob = 1;
                        % record change magnitude
                        rec_cg(num_fc).magnitude = median(v_dif_mag,1);
                        
                        % identified and move on for the next functional curve
                        num_fc = num_fc + 1;
                        % start from i+1 for the next functional curve
                        i_start = i + 1;
                        % start training again
                        BL_train = 0;
                        
                    elseif vec_mag(1) > Tmax_cg % false change
                        % remove noise
                        clrx(i+1) = [];
                        clry(i+1,:) = [];
                        i=i-1; % stay & check again after noise removal
                    end
                end % end of continuous monitoring
            end % end of checking basic requrirements
            
            % move forward to the i+1th clear observation
            i=i+1;
        end % end of while iterative
        
        % Two ways for processing the end of the time series
        if BL_train == 1
            % 1) if no break find at the end of the time series
            % define probability of change based on conse
            for i_conse = conse:-1:1
                if vec_mag(i_conse) <= T_cg
                    % the last stable id
                    id_last = i_conse;
                    break;
                end
            end
            
            % update change probability
            rec_cg(num_fc).change_prob = (conse-id_last)/conse;
            % update end time of the curve
            rec_cg(num_fc).t_end=clrx(end-conse+id_last);
            
            if conse > id_last % > 1
                % update time of the probable change
                rec_cg(num_fc).t_break = clrx(end-conse+id_last+1);
                % update magnitude of change
                rec_cg(num_fc).magnitude = median(v_dif_mag(id_last+1:conse,:),1);
            end
            
        elseif BL_train == 0
            % 2) if break find close to the end of the time series
            % Use [conse,min_num_c*n_times+conse) to fit curve
            
            if num_fc == rec_fc
                % first curve
                i_start = 1;
            else
                i_start = find(clrx >= rec_cg(num_fc-1).t_break);
                i_start = i_start(1);
            end
            
            % Tmask
            if length(clrx(i_start:end)) > conse
                blIDs = autoTmask(clrx(i_start:end),clry(i_start:end,[num_B1,num_B2]),...
                    (clrx(end)-clrx(i_start))/num_yrs,adj_rmse(num_B1),adj_rmse(num_B2),T_const);
                
                % update i_span after noise removal
                i_span = sum(~blIDs);
                
                IDs = i_start:length(clrx); % all IDs
                rmIDs = IDs(blIDs(1:end-conse) == 1); % IDs to be removed
                
                % remove noise pixels between i_start & i
                clrx(rmIDs) = [];
                clry(rmIDs,:) = [];
            end
            
            % enough data
            if length(clrx(i_start:end)) >= conse
                % defining computed variables
                fit_cft = zeros(max_num_c,nbands-1);
                % rmse for each band
                rmse = zeros(nbands-1,1);
                % end of fit qa = 20
                qa = 20;
                
                for i_B = 1:nbands-1
                    [fit_cft(:,i_B),rmse(i_B)] = ...
                        autoTSFit(clrx(i_start:end),clry(i_start:end,i_B),min_num_c);
                end
                
                % record time of curve start
                rec_cg(num_fc).t_start = clrx(i_start);
                % record time of curve end
                rec_cg(num_fc).t_end=clrx(end);
                % record break time
                rec_cg(num_fc).t_break = 0;
                % record postion of the pixel
                rec_cg(num_fc).pos = (nrows-1)*ncols+i_ids;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = length(clrx(i_start:end));
                % record fit category
                rec_cg(num_fc).category = qa + min_num_c;
                % record change magnitude
                rec_cg(num_fc).magnitude = zeros(1,nbands-1);
            end
        end
    end % end of if sum(idgood) statement
end % end of for i_ids loop
save(fullfile(dir_l,n_rst,['record_change',num2str(nrows)]),'rec_cg');
end % end of function
