function autoShowLandDistMap(varargin)
% This function is used to provde 1) disturbance maps for each year and 2)
% most recent disturbance

% Specific parameters
% ------------------------
%   'CCDCDir'     Directory of input data.  Default is the path to
%                        the current folder.
%   'StartYear'          Start year.
%   'EndYear'            End year.
%
% autoShowLandDistMap('StartYear', 2010,'EndYear', 2015)
% will uotput the disturbance maps between 2010 and 2015 as well as the
% most recent disturbance (within all observations)
%
% Version 1.00 No disturbance class in the cover map (03/29/2018)

 %% get parameters from inputs
% where the all Landsat zipped files are
dir_cur = pwd;
p = inputParser;
p.FunctionName = 'paras';

addParameter(p,'CCDCDir',dir_cur);
addParameter(p,'StartYear',0);
addParameter(p,'EndYear',0);
 % request user's input
parse(p,varargin{:});
dir_cur=p.Results.CCDCDir;
start_year=p.Results.StartYear;
end_year=p.Results.EndYear;

% version of CCDC change
try
    vr = varead(fullfile(dir_cur,'CCDC_Change_log.txt'),'Version');
catch
    vr= '0'; % indicatting unknown version.
end

% get image parameters automatically`
imf=dir(fullfile(dir_cur,'L*')); % folder names
[nrows,ncols,nbands,jiUL,res,zc,~] = autoPara(imf);
l_dir = dir_cur;

% INPUTS:
all_yrs = start_year:end_year;%1985:2015; % all of years for producing maps

% dimension and projection of the image
jiDim = [ncols,nrows];
% max number of maps
max_n = length(all_yrs);
% slope threshold
t_min = 0; % 0.01 change in surf ref in ten years

% produce disturbance map
LandDistMap = 9999*ones(nrows,ncols,max_n,'uint16'); % disturbance magnitude
% most recent accumulated map
AccuYearMap = 9999*ones(nrows,ncols,'uint16'); % disturbance magnitude

% make Predict folder for storing predict images
n_map = 'CCDCMap';
if isempty(dir(fullfile(dir_cur,n_map)))
    mkdir(n_map);
end

% cd to the folder for storing recored structure
% cd(v_input.name_rst);
n_str = 'TSFitMap';
imf = dir(fullfile(dir_cur,n_str,'record_change*')); % folder names
num_line = size(imf,1);

for line = 1:num_line
    
    % show processing status
    if line/num_line < 1
        fprintf('Processing %.2f percent\r',100*(line/num_line));
    else
        fprintf('Processing %.2f percent\n',100*(line/num_line));
    end
    
    % load one line of time series models
    load(fullfile(dir_cur,n_str,imf(line).name)); %#ok<LOAD>
    
    % postions
    pos = [rec_cg.pos];
    
    % continue if there is no model available
    l_pos = length(pos);
    if l_pos == 0
        continue
    end
    
    % break time
    t_break = [rec_cg.t_break];
    % change probability
    change_prob = [rec_cg.change_prob];
    % change vector magnitude
    mag = [rec_cg.magnitude];
    % reshape magnitude
    mag = reshape(mag,nbands-1,[]);
    
    
    for i = 1:l_pos
        % get row and col
        [I,J] = ind2sub(jiDim,pos(i));
        
        % initialize pixels have at least one model
        if sum(LandDistMap(J,I,:) == 9999) == max_n
            % write doy to ChangeMap
            LandDistMap(J,I,:) = 0;
            AccuYearMap(J,I,:) = 0;
        end
        
        if change_prob(i) == 1
            [break_type,break_year,break_doy] = label_dist_vec(t_break(i),t_min,mag(:,i));
            
            if break_type == 3 % abrupt distrubance
                % get the band number for abrupt disturbance
                n_band = all_yrs == break_year;
                % magnitude of disturbance
                LandDistMap(J,I,n_band) = break_doy;
                % update accumulated by year map
                AccuYearMap(J,I) = break_year;
            end
        end
    end
end

% Scence or ARD
if zc >=1 && zc <= 60
    % write ENVI files
    enviwrite_bands(fullfile(l_dir,n_map,['LandDistMap',char(vr)]),LandDistMap,'uint16',res,jiUL,'bsq',zc,all_yrs);
    enviwrite_bands(fullfile(l_dir,n_map,['AccuYearMap',char(vr)]),AccuYearMap,'uint16',res,jiUL,'bsq',zc,start_year*100000+end_year);
else
    % confirmed change
    ARD_enviwrite_bands(fullfile(l_dir,n_map,['LandDistMap',char(vr)]),LandDistMap,'uint16','bsq',all_yrs);
    ARD_enviwrite_bands(fullfile(l_dir,n_map,['AccuYearMap',char(vr)]),AccuYearMap,'uint16','bsq',start_year*100000+end_year);
end