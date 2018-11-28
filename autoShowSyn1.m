function autoShowSyn1(varargin)
% This function is used to provde synthetic data for a certain date (mm/dd/yyyy).
% two digits for QA band
%
% Specific parameters
% ------------------------
%   'InDir'     Directory of CCDC change detection results.  
%                        Default is the path to current folder.
%   'OutDir'    Directory of the model output.  Default is the 
%                        same as InputDirectory.
%   'Year'               Year.
%   'Month'              Month.
%   'Day'                Day.
%   'DOY'                Day of year.

% Examples:
% ------------------------
% autoShowSyn1('Year',2010,'Month',7,'Day',1)
% will give the synthetic data at July 1,2010.
%
% autoShowSyn1('Year',2010,'DOY',200)
% will give the synthetic data at 200th day,2010.
%

 %% get parameters from inputs
% where the all CCDC change detection results are
dir_cur = pwd;
% where the output files are
dir_out = '';
p = inputParser;
p.FunctionName = 'paras';

addParameter(p,'InDir',dir_cur);
addParameter(p,'OutDir',dir_out);
addParameter(p,'Year',0);
addParameter(p,'Month',0);
addParameter(p,'Day',0);
addParameter(p,'DOY',0);

parse(p,varargin{:});
dir_cur = p.Results.InDir;
dir_out = p.Results.OutDir;

yy=p.Results.Year;
mm=p.Results.Month;
dd=p.Results.Day;
doy=p.Results.DOY;
if isempty(dir_out)
    dir_out = dir_cur;
end

if yy>0
    if mm>0&&dd>0 % using yyyy/mm/dd
        % converted to julian date
        j_date=datenum(yy,mm,dd);
        j0_date=datenum(yy,1,0);
        doy = j_date-j0_date;
        % date for show i.e. 19990815
        %s_date=yy*10000+mm*100+dd;
        s_date = yy*1000 + doy;
    elseif doy>0 % using yyyy DOY
        j_date = datenum(yy,1,0) + doy;
        s_date = yy*1000 + doy;
    end
end
% v_input = ccdc_Input(dir_cur);

% 1. inputs: the time for predicted synthetic data
% yy=2000;
% mm=12;
% dd = 1;
% OR
% yy=2000;
% doy=120;

% nbands = v_input.nbands - 1;% number of bands in the image
% ncoefs = v_input.num_c;% number of coefficients
% 
% % dimension and projection of the image
% % cd(v_input.l_dir);
% nrows = v_input.ijdim(1);
% ncols = v_input.ijdim(2);
% jiDim = [ncols,nrows];
% jiUL = v_input.jiul;
% res = v_input.resolu;
% zc = v_input.zc;

% get num of total folders start with "L"
imf=dir(fullfile(dir_cur,'L*')); % folder names
% filter for Landsat folders
imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
imf = [imf{:}];
imf = vertcat(imf{:});
% sort according to yeardoy
yeardoy = str2num(imf(:, 10:16));
[~, sort_order] = sort(yeardoy);
imf = imf(sort_order, :);
% name of the first stacked image
filename = dir(fullfile(dir_cur,imf(1,:),'L*MTLstack')); 
% read in dimension and zone number of the data
[jiDim,jiul,resolu,zc] = envihdrread(fullfile(dir_cur,imf(1,:),filename.name));

n_syn = ['LS',imf(1,4:9),num2str(s_date)];

nbands = 8; %original band amounts
nbands = nbands - 1; % show change for all the bands except fmask

ncoefs = 8; %max num of coefficients
nrows = jiDim(2);
ncols = jiDim(1);

% produce synthetic data (1, 2, 3, 4, 5, 7, and QA)
SR = -9999*ones(nrows,ncols,nbands,'int16'); % 1. change here for excluding thermal

% make Predict folder for storing predict images
%n_pre = v_input.name_pre;%'PredictAll';
n_pre = 'PredictAll';

if isempty(dir(fullfile(dir_cur,n_pre)))
    mkdir(fullfile(dir_cur,n_pre));
end

% folder saved everything
%n_rst = v_input.name_rst;
n_rst = 'TSFitMap';
%cd(n_rst);% TSFitMap

imf=dir(fullfile(dir_cur,n_rst,'record_change*')); % folder names
num_line = size(imf,1);

%%
for line = 1:num_line
    %fprintf('Processing %.2f percent\r',100*(line/num_line));
    load(fullfile(dir_cur,n_rst,imf(line).name));
    
    % postions & coefficients
    pos = [rec_cg.pos];
    % continue if there is no model available
    l_pos = length(pos);
    if l_pos == 0
        continue;
    end
    % get coefficients matrix
    coefs = [rec_cg.coefs];
    % reshape coefs
    coefs = reshape(coefs,nbands*ncoefs,l_pos);
    % get category matrix
    category = [rec_cg.category];
    % take the first digit (number of coefficients)
    category = category - 10*floor(category/10);
    
    % model start, end, & break time for prediction
    model_start = [rec_cg.t_start];
    model_end = [rec_cg.t_end];
    model_break = [rec_cg.t_break];
    % model on the right
    ids_right = model_start > j_date;
    % model on the left
    ids_left = (model_end < j_date & model_break == 0)|(model_break <= j_date & model_break > 0);
    % id within model interval
    ids_middle = model_start <= j_date & ((model_end >= j_date & model_break == 0) | (model_break > j_date & model_break > 0));

    % position for model in the middle
    pos_middle = pos(ids_middle);
    % coefficients for model in the middle
    coefs_middle = coefs(:,ids_middle);
    % category for model in the middle
    category_middle = category(ids_middle);
    
    % positions for the nearest model on the right
    pos_right = pos(ids_right);
    [pos_near_right,ids_near_right] = unique(pos_right,'first');
    % coefficients for the nearest model on the right
    coefs_right = coefs(:,ids_right);
    coefs_near_right = coefs_right(:,ids_near_right);
    % category for the nearest model on the right
    category_right = category(ids_right);
    category_near_right = category_right(ids_near_right);
    
    % postions for the nearest model on the left
    pos_left = pos(ids_left);
    [pos_near_left,ids_near_left] = unique(pos_left,'last');
    % coefficients for the nearest model on the left
    coefs_left = coefs(:,ids_left);
    coefs_near_left = coefs_left(:,ids_near_left);
    % category for the nearest model on the left
    category_left = category(ids_left);
    category_near_left = category_left(ids_near_left);
    
    % pass if there is no nearest model on the left 
    l_pos=length(pos_near_left);
    if l_pos > 0  
        % providing predictions
        for i=1:l_pos
            [I,J]=ind2sub(jiDim,pos_near_left(i));
            for j_b=1:nbands-1 % excluding thermal band
                SR(J,I,j_b)=autoTSPred(j_date,coefs_near_left(((j_b-1)*ncoefs+1):j_b*ncoefs,i));
            end
            % QA band
            SR(J,I,end) = 20 + category_near_left(i); % model forward predicted values         
        end
    end
    
    % pass if there is no nearest model on the right 
    l_pos=length(pos_near_right);
    if l_pos > 0  
        % providing predictions
        for i=1:l_pos
            [I,J]=ind2sub(jiDim,pos_near_right(i));
            for j_b=1:nbands-1 % excluding thermal band
                SR(J,I,j_b)=autoTSPred(j_date,coefs_near_right(((j_b-1)*ncoefs+1):j_b*ncoefs,i));
            end
            % QA band
            SR(J,I,end) = 10 + category_near_right(i); % model backward predicted values         
        end
    end
    
    % pass if there is no nearest model in the middle 
    l_pos=length(pos_middle);
    if l_pos > 0  
        % providing estimations
        for i=1:l_pos
            [I,J]=ind2sub(jiDim,pos_middle(i));
            for j_b=1:nbands-1 % excluding thermal band
                SR(J,I,j_b)=autoTSPred(j_date,coefs_middle(((j_b-1)*ncoefs+1):j_b*ncoefs,i));
            end
            % QA band
            SR(J,I,end) = 0 + category_middle(i); % model estimated values         
        end
    end    
end

n_bands = {'Band 1','Band 2','Band 3','Band 4','Band 5','Band 7','QA'};
% n_bands = 1:7;
% write ENVI files
ARD_enviwrite_bands_n(fullfile(dir_out,n_pre,n_syn),SR,'int16','bip',n_bands,dir_cur,dir_out);
end





