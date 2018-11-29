function autoShowLandDistMap(varargin)
% This function is used to provde disturbance maps for each year
% Version 1.00 No disturbance class in the cover map (03/29/2018)
% vr = varead('COLD_log.txt','Version');
%
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


% get image parameters automatically
% get parameters from inputs
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

imf = dir(fullfile(dir_cur,'L*')); % folder names

% filter for Landsat folders
imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
imf = [imf{:}];
imf = vertcat(imf{:});
% name of the first stacked image
filename = dir(fullfile(dir_cur,imf(1,:),'L*stack'));
% read in ENVI hdr
info = read_envihdr(fullfile(dir_cur,imf(1,:),[filename.name,'.hdr']));
% provide values from info
nrows = info.lines;
ncols = info.samples;
nbands = info.bands;
% get current directory
l_dir = dir_cur;

% INPUTS:
all_yrs = start_year:end_year;% all of years for producing maps
% dimension and projection of the image
jiDim = [ncols,nrows];
% max number of maps
max_n = length(all_yrs);
% slope threshold
t_min = -200; % 0.02 change in surf ref 

% produce disturbance map
LandDistMap = 9999*ones(nrows,ncols,max_n,'uint16'); % disturbance magnitude
% most recent accumulated map
AccuYearMap = 9999*ones(nrows,ncols,'uint16'); % disturbance magnitude
% most recent accumulated map
AccuTypeMap = 9999*ones(nrows,ncols,'uint8'); % disturbance magnitude

% make Predict folder for storing predict images
n_map = 'CCDCMap';
if isempty(dir(fullfile(dir_cur,n_map)))
    mkdir(fullfile(dir_cur,n_map));
end

% cd to the folder for storing recored structure
% cd(v_input.name_rst);
n_str = 'TSFitMap';
imf = dir(fullfile(dir_cur,n_str,'record_change*')); % folder names
num_line = size(imf,1);

for line = 1300:num_line
    
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
    % coefficients
    coefs = [rec_cg.coefs];
    coefs = reshape(coefs,8,nbands-1,[]);
    
    
    for i = 1:l_pos - 1
        % get row and col
        [I,J] = ind2sub(jiDim,pos(i));
        
        % initialize pixels have at least one model
        if sum(LandDistMap(J,I,:) == 9999) == max_n
            % write doy to ChangeMap
            LandDistMap(J,I,:) = 0;
            AccuYearMap(J,I,:) = 0;
        end
        
        if change_prob(i) == 1
            [break_type,break_year,break_doy] = label_dist_type(coefs(:,:,i),t_break(i),t_min,mag(:,i),coefs(:,:,i+1));
            % [break_type,break_year,break_doy] = label_dist_type(tst(j).coefs,tst(j).t_break,-200,tst(j).magnitude,tst(j+1).coefs);
            
            if break_type > 1 % land distrubance
                % get the band number for abrupt disturbance
                n_band = all_yrs == break_year;
                % magnitude of disturbance
                LandDistMap(J,I,n_band) = break_doy;
                % update accumulated by year map
                AccuYearMap(J,I) = break_year;
            end
            
            % update accumulated type map
            AccuTypeMap(J,I) = break_type;
        end
    end
end

rs_imwrite_bands(LandDistMap,fullfile(l_dir,n_map,'LandDistMap'),info,all_yrs); 
rs_imwrite_bands(AccuYearMap,fullfile(l_dir,n_map,'AccuYearMap'),info,start_year*100000+end_year); 
rs_imwrite_bands(AccuTypeMap,fullfile(l_dir,n_map,'AccuTypeMap'),info,'Disturb Type');

% % Scence or ARD
% if zc >=1 && zc <= 60
%     % write ENVI files
%     enviwrite_bands([l_dir,'/',n_map,'/LandDistMap',char(vr)],LandDistMap,'uint16',res,jiUL,'bsq',zc,all_yrs);
%     enviwrite_bands([l_dir,'/',n_map,'/AccuYearMap',char(vr)],AccuYearMap,'uint16',res,jiUL,'bsq',start_year*100000+end_year);
% else
%     % confirmed change
%     ARD_enviwrite_bands([l_dir,'/',n_map,'/LandDistMap',char(vr)],LandDistMap,'uint16','bsq',all_yrs);
%     ARD_enviwrite_bands([l_dir,'/',n_map,'/AccuYearMap',char(vr)],AccuYearMap,'uint16','bsq',start_year*100000+end_year);
% end