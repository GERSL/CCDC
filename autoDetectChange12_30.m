function autoDetectChange12_30(task,ntasks,varargin)
% Matlab code for Continous Change Detection (for standalone version)

% Specific parameters
% ------------------------
%   'CCDCDir'     Directory of CCDC change detection results.  
%                 Default is the path to current folder.
% Inputs:
% T_cg:   change probability threshold
% conse:  number of consecutive observations
% max_c:  maximum number of coefficients used
% task:   The task CPU used
% ntasks: number of CPUs used
% profile on;
% Get variable and path 
dir_l = pwd; %'/projectnb/landsat/projects/LCMS/4530/images/'; 
p = inputParser;
p.FunctionName = 'paras';
addParameter(p,'CCDCDir',dir_l);

parse(p,varargin{:});
dir_l = p.Results.CCDCDir;

%% Tuning variables for continuous change detection
% probability for detecting surface change
T_cg = 0.99; % (default) 
% number of consecutive observation 
conse = 6; % (default)
% number of maximum coefficients
max_c = 8; % (default)

%% Constants: 
% Bands for detection change
B_detect = 2:6;
% Tmasking of noise
Tmax_cg = chi2inv(1-1e-6,length(B_detect)); 
% version of CCDC
ccdc_v = 12.30;

% read in image dimension from reference image
% imf = dir('L*');
% filename = [dir_l,'/example_img'];

% get num of total folders start with "L"
imf=dir(fullfile(dir_l,'L*')); % folder names
% number of images
num_imgs = size(imf,1);
% filter for Landsat folders
imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
imf = [imf{:}];
imf = vertcat(imf{:});
% sort according to yeardoy
yeardoy = str2num(imf(:, 10:16));
[~, sort_order] = sort(yeardoy);
imf = imf(sort_order, :);
% name of the first stacked image
filename = dir(fullfile(dir_l,imf(1,:),'L*MTLstack')); 
% read in dimension and zone number of the data
[jiDim,~,~,zc] = envihdrread(fullfile(dir_l,imf(1,:),filename.name));
% dimension of image [row,col]
ijdim = [jiDim(2),jiDim(1)];
% number of nrows processed
nrows = ijdim(1);
% number of pixels procesed per line
ncols = ijdim(2);
% total ouput bands (1-5,7,6,cfmask)
nbands = 8;

% folder name of all CCDC results 
name_rst = 'TSFitMap';
% make TSFitMap folder for storing coefficients
if isempty(dir(fullfile(dir_l,name_rst)))
    mkdir(fullfile(dir_l,name_rst));
end

% log for CCDC Change paramters and versions
% report only for the first task
if task == 1
    fileID = fopen(fullfile(dir_l,'CCDC_Change_log.txt'),'w');
    % write location of image stack
    fprintf(fileID,'Image location = %s\r\n',dir_l);
    % write number of images used
    fprintf(fileID,'Number of images = %d\r\n',num_imgs);
    % change probability threshold
    fprintf(fileID,'Change Probability Threshold = %.2f\r\n',T_cg);
    % # of consecutive observations
    fprintf(fileID,'Number of Consecutive Observations = %d\r\n',conse);
    % Maximum number of coefficients
    fprintf(fileID,'Maximum Number of Coefficients = %d\r\n',max_c);
    % CCDC Version
    % fprintf(fileID,'CCDC Change for ARD Version = %.2f\r\n',ccdc_v);
    fprintf(fileID,'CCDC Change for Single Scene Version = %.2f\r\n',ccdc_v);
    % updates
    fprintf(fileID,'******************************************************************************************************\r\n');
    fprintf(fileID,'Revisions: $ Date: 11/20/2015 $ Copyright: Zhe Zhu\r\n');
    fprintf(fileID,'Version 12.30  Fixed a bug for pixels without minimum observations (11/20/2015)\r\n');
    fprintf(fileID,'Version 12.29  Modify fit for perennial snow and Fmask failed pixels (09/20/2015)\r\n');
    fprintf(fileID,'Version 12.28  Fix a bug for missing values in land cover maps (09/16/2015)\r\n');
    fprintf(fileID,'Version 12.27  Fix a bug for persistent snow and falied Fmask pixels (06/17/2015)\r\n');
    fprintf(fileID,'Version 12.26  Connect time for all models (05/28/2015)\r\n');
    fprintf(fileID,'Version 12.25  Fix a bug in snow percent (05/19/2015)\r\n');
    fprintf(fileID,'Version 12.24  Change T_const in Tmask (03/31/2015)\r\n');
    fprintf(fileID,'Version 12.23  Update iteratively before 24 observations (03/22/2015)\r\n');
    fprintf(fileID,'Version 12.22  Adjust mini RMSE based on temporal variability (03/22/2015)\r\n');
    fprintf(fileID,'Version 12.21  Add more categories and update i_start in the end (03/14/2015)\r\n');
    fprintf(fileID,'Version 12.20  Convert BT to from K to C before analysis (03/12/2015)\r\n');
    fprintf(fileID,'Version 12.19  Fit for permanent snow if it is > 75 percent (03/12/2015)\r\n');
    fprintf(fileID,'Version 12.18  No change detection if clear observation < 25 percent (03/12/2015)\r\n');
    fprintf(fileID,'Version 12.17  Use median value for very simple model & change magnitude (02/24/2015)\r\n');
    fprintf(fileID,'Version 12.16  Finding changes in all water pixels (02/24/2015)\r\n'); 
    fprintf(fileID,'Version 12.15  Use the original multitemporal cloud mask (02/15/2015)\r\n');
    fprintf(fileID,'Version 12.14  Do not need example_img in images folder (02/09/2015)\r\n');
    fprintf(fileID,'Version 12.13  More infromation in "category" (11/10/2014)\r\n');    
    fprintf(fileID,'Version 12.12  Fit for pixels where Fmask fails (11/09/2014)\r\n');
    fprintf(fileID,'Version 12.11  Fix a bug in num_fc (11/09/2014)\r\n');
    fprintf(fileID,'Version 12.10  Better multietmporal cloud detection at the beginning (11/06/2014)\r\n');
    fprintf(fileID,'Version 12.09  Detect change for land pixel (water/snow speical case) (10/31/2014)\r\n');
    fprintf(fileID,'Version 12.08  Speed up by reducing time for RMSE and model computing (10/17/2014)\r\n');
    fprintf(fileID,'Version 12.07  mini rmse should be larger than 10 percent of the mean (10/13/2014)\r\n');
    fprintf(fileID,'Version 12.06  Fit model again when there are a 33 percent more in time (10/08/2014)\r\n');
    fprintf(fileID,'Version 12.05  Use subset of bands (2-6) for detecting surface change (10/01/2014)\r\n');
    fprintf(fileID,'Version 12.04  Only apply Tmask during model initialization (09/29/2014)\r\n');
    fprintf(fileID,'Version 12.03  Use subset of bands (3-5) for detecting surface change (09/01/2014)\r\n');
    fprintf(fileID,'Version 12.02  Fix a bug in model intialization (08/14/2014)\r\n');
    fprintf(fileID,'Version 12.01  Use subset of bands (3-6) to reduce atmosphere influence (08/04/2014)\r\n');
    fprintf(fileID,'Version 12.00  Detect change based on probability (07/19/2014)\r\n');
    fprintf(fileID,'Version 11.06  No need to change folder name & faster in speed (06/06/2014)\r\n');
    fprintf(fileID,'Version 11.05  Improve calculation of temporally adjusted RMSE (04/23/2014)\r\n');
    fprintf(fileID,'Version 11.04  Revise "rec_cg.category" for different fit processes (04/01/2014)\r\n');
    fprintf(fileID,'Version 11.03  Add "rec_cg.magnitude" as change magnitude indicator (04/01/2014)\r\n');
    fprintf(fileID,'Version 11.02  Change very simple fit with mean value (04/01/2014)\r\n');
    fprintf(fileID,'Version 11.01  Do not need metadata in the image folder to run CCDC (03/25/2014)\r\n');
    fprintf(fileID,'Version 11.00  Use change vector magnitude as hreshold for change (03/25/2014)\r\n');
    fprintf(fileID,'Version 10.13  Use subset of bands (1-6) to reduce atmosphere influence (01/31/2014)\r\n');
    fprintf(fileID,'Version 10.12  More accurate number of days per year "num_yrs" (01/30/2014)\r\n');
    fprintf(fileID,'Version 10.11  RMSE updates with time series fit (01/26/2014)\r\n');
    fprintf(fileID,'Version 10.10  Update temperature extreme in recent studies (01/16/2014)\r\n');
    fprintf(fileID,'Version 10.09  Find break in max value in any of the band (01/08/2014)\r\n');
    fprintf(fileID,'Version 10.08  Add very simple fit (median) for start & end of timeseries (10/21/2013)\r\n');
    fprintf(fileID,'Version 10.07  Better multitemporal cloud detection (10/19/2013)\r\n');
    fprintf(fileID,'Version 10.06  Add "Tmax_cg" for last step noise removal (10/18/2013)\r\n');
    fprintf(fileID,'Version 10.05  Use subset of bands (2-6) to avoid atmosphere influences (10/18/2013)\r\n');
    fprintf(fileID,'Version 10.04  Let dynamic fitting for pixels at the beginning (09/23/2013)\r\n');
    fprintf(fileID,'Version 10.03  Able to detect change at the verying beginning (09/06/2013)\r\n');
    fprintf(fileID,'Version 10.02  Add mini years "mini_yrs" in model intialization (09/03/2013)\r\n');
    fprintf(fileID,'Version 10.01  Reduce time for calcuating "v_dif" (09/02/2013)\r\n');
    fprintf(fileID,'Version 10.00  Fit for beginning and end of the time series (08/31/2013)\r\n');
    fprintf(fileID,'Version 9.09   Only fit more than 50 percent of Landat images overlap area (08/28/2013)\r\n');
    fprintf(fileID,'Version 9.08   Force model fit for persistent snow pixels (08/27/2013)\r\n');
    fprintf(fileID,'Version 9.07   Add "rec_cg.category", "rec_cg.change_prob", and "recc_cg.num_obs" (08/20/2013)\r\n');
    fprintf(fileID,'Version 9.06   Remove mininum rmse "mini" and minimum years "mini_yrs" (08/16/2013)\r\n');
    fprintf(fileID,'Version 9.05   Model gets more coefficients with more observations (08/16/2013)\r\n');
    fprintf(fileID,'Version 9.04   Fix a bug in calculating temporally adjusted rmse (08/01/2013)\r\n');
    fprintf(fileID,'Version 9.03   Fit curve again after one year (03/28/2013)\r\n');
    fprintf(fileID,'Version 9.02   Use "mini = T_const/T_cg" for small rmse cases (03/26/2013)\r\n');
    fprintf(fileID,'Version 9.01   Remove out of range pixels before time series analysis (02/09/2013)\r\n');
    fprintf(fileID,'Version 9.00   Using 8 coefficients and lasso fit (02/01/2013)\r\n');
    fprintf(fileID,'Version 8.04   Use "max v_slope" instead of "average v_slope" (01/16/2013)\r\n');
    fprintf(fileID,'Version 8.03   Start initialization when "time_span" > 1 year (01/16/2013)\r\n');
    fprintf(fileID,'Version 8.02   Fix a bug in not fitting models at the begining (01/16/2013)" (08/20/2013)\r\n');
    fprintf(fileID,'Version 8.01   Fix a bug in counting "i" and "i_span"(01/13/2013)\r\n');
    fprintf(fileID,'Version 8.00   Temporally changing RMSE (01/09/2013)\r\n');
    fprintf(fileID,'Version 7.03   Continuous Change Detection and Classification (CCDC) (07/11/2012)\r\n');
    fprintf(fileID,'Version 1.00   Continous Monitoring of Forest Disturbance Algorithm (CMFDA) (07/13/2010) \r\n');
    fprintf(fileID,'******************************************************************************************************\r\n');
    fclose(fileID);
end

% adjust threshold based on chi-squared distribution
% T_cg = chi2inv(T_cg,1)^0.5; 
T_cg = chi2inv(T_cg,length(B_detect));

% prepare the irows for idn_cpu for ALL rows
irows = zeros(1,1);
i = 0;
while task + ntasks*i <= nrows % process all lines
   irows(i+1) = task + ntasks*i;
   i = i+1;
end

% % only calculate a few lines
% % start with line_start
% line_start = 3700;
% line_end = 4000;
% 
% irows(irows < line_start) = [];
% irows(irows > line_end) = [];

for i = 1:length(irows)
    % Check whether record_change already exist for row == i
    try
        % if successfully loaded => skip
        load(fullfile(dir_l,name_rst,['record_change',sprintf('%d',irows(i)),'.mat']));
    catch me
        % not exist or corrupt
        fprintf('Processing the %dth row\n',irows(i));        
        % Continous Change Detection Done for a line of timeseries pixels
        % TrendSeasonalFit_v12_30ARDLine(dir_l,name_rst,ncols,irows(i),T_cg,Tmax_cg,conse,max_c,nbands,B_detect);
        TrendSeasonalFit_v12_30Line(dir_l,name_rst,ncols,irows(i),T_cg,Tmax_cg,conse,max_c,nbands,B_detect);
    end
end

% profile viewer
exit
end % end of function 