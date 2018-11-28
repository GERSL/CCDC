function autoShowClassMap(varargin)
% This function is used to provde all land classification maps for each year
% Results for LCMAP
% Version 1.2  Acquire map from annual classfication map (04/09/2015)
% Version 1.1  Add disturbance in the annual map (04/01/2015);
% Version 1.0  No disturbance class in the cover map (11/06/2015)
% Tools
%
%
% Specific parameters
% ------------------------
%   'CCDCDir'     Directory of input data.  Default is the path to
%                        the current folder.
%   'StartYear'          Start year.
%   'EndYear'            End year.
%   'Year'               Year. Only a certain date's map will be outputted
%                        when inputing this.
%   'Month'              7
%   'Day'                1
%
%
% Examples:
% ------------------------
% autoShowClassMap('StartYear',2010,'EndYear',2015)
% will give the classifiaction maps at July 1 between 2010 and 2015.
%
% autoShowClassMap('StartYear',2000,'EndYear',2015,'Month',12,'Day',1)
% will give the classifiaction maps at December 1 between 2010 and 2015.
%
% autoShowClassMap('Year',2000,'Month',7,'Day',1)
% will give the classifiaction map at July 1, 2000.
%
   %% get parameters from inputs
    % where the all Landsat zipped files are
    dir_cur = pwd;
    p = inputParser;
    p.FunctionName = 'paras';
    
    addParameter(p,'CCDCDir',dir_cur);
    
    % default date is 1th, July.
    addParameter(p,'StartYear',0);
    addParameter(p,'EndYear',0);
    addParameter(p,'Year',0);
    addParameter(p,'Month',7);
    addParameter(p,'Day',1);
    
    % request user's input
    parse(p,varargin{:});
    dir_cur=p.Results.CCDCDir;
    
    start_year=p.Results.StartYear;
    end_year=p.Results.EndYear;
    oneyear=p.Results.Year;
    
    mm=p.Results.Month;
    dd=p.Results.Day;
    
    % give certain year
    if oneyear>0
        start_year = oneyear;
        end_year = oneyear;
    end
    
    % optional
%     addpath('~/ccdc');
    
    cd(dir_cur);

    % INPUTS:
    all_yrs = start_year:end_year; % all of years for producing maps
    v_input = ccdc_Inputs;
    
    % dimension and projection of the image
    nrows = v_input.ijdim(1);
    ncols = v_input.ijdim(2);
    jiDim = [ncols,nrows];
    jiUL = v_input.jiul;
    res = v_input.resolu;
    zc = v_input.zc;
    % number of coefficients
    num_c = v_input.num_c;
    % number of bands
    nbands = v_input.nbands;
    % max number of maps
    max_n = length(all_yrs);
    % all julidan dates
%     mm = 7;
%     dd = 1;
    jul_d = datenummx(all_yrs,mm,dd);
    jul_start = datenummx(all_yrs,1,1);
    jul_end = datenummx(all_yrs,12,31);

    % produce land cover map
    CoverMap = 255*ones(nrows,ncols,max_n,'uint8'); % Trends categories (0~11)
    % produce land cover map QA (unsupervised emsemble margin)
    CoverQAMap = 255*ones(nrows,ncols,max_n,'uint8'); % Trends categories (0~100)

    % make Predict folder for storing predict images
    n_map=v_input.name_map;% 'CCDCMap';
    if isempty(dir(n_map))
        mkdir(n_map);
    end

    % cd to the folder for storing recored structure
    cd(v_input.name_rst);

    imf = dir(fullfile(dir_cur,'record_change*')); % folder names
    num_line = size(imf,1);
    line_pt = 0;
    %%
    for line = 1:num_line
        if 100*(line/num_line) - line_pt > 1
            fprintf('Processing %.0f percent\n',ceil(100*(line/num_line)));
            line_pt = 100*(line/num_line);
        end

        % load one line of time series models
        load(fullfile(dir_cur,imf(line).name));

        % postions
        pos = [rec_cg.pos];

        % continue if there is no model available
        l_pos = length(pos);
        if l_pos == 0
            continue
        end

        % matrix of each component
        % start time
        t_start = [rec_cg.t_start];
        % end time
        t_end = [rec_cg.t_end];
        % break time
        t_break = [rec_cg.t_break];
        % class
        class = [rec_cg.class];
        % classification QA
        class_qa = [rec_cg.classQA];

% %         % put disturbed class to grass/shrub
% %         class(class == 10) = 7;

        % For Condition & Cover Map
        for i = 1:l_pos
            % get row and col
            [I,J] = ind2sub(jiDim,pos(i));

            % initialize pixels have at least one model
            if sum(CoverMap(J,I,:) == 255) == max_n
                % write land cover to CoverMap
                CoverMap(J,I,:) = 0;
                % write land cover to CoverQAMap
                CoverQAMap(J,I,:) = 0;
            end

            % give disturbed class for each year        
            % year (band) the curve belongs to
            n_band = jul_d >= t_start(i) & (jul_d <= t_end(i) | jul_d < t_break(i));        
            % write land cover to CoverMap
            CoverMap(J,I,n_band) = class(i);
            % write land cover to CoverQAMap
            CoverQAMap(J,I,n_band) = class_qa(i);        

            % give next land cover category for gaps
            if i > 1
                if pos(i) == pos(i-1) % same location
                    n_dist = jul_d < t_start(i) & jul_d >= t_break(i-1);
                    % find the first nonzero value from class(:,i);
                    idn = find(class(:,i)>0);
                    if ~isempty(idn)
                        % write land cover (distubed) to CoverMap
                        CoverMap(J,I,n_dist) = class(idn(1),i);
                    end
                end
            end
        end
    end
    cd ..

    enviwrite_bands(fullfile(v_input.l_dir,n_map,'CoverMap'),CoverMap,'uint8',res,jiUL,'bsq',zc,all_yrs);
    clear CoverMap;
    enviwrite_bands(fullfile(v_input.l_dir,n_map,'CoverQAMap'),CoverQAMap,'uint8',res,jiUL,'bsq',zc,all_yrs);
    clear CoverQAMap;

% %     % write ENVI files for ARD
% %     % Cover Map
% %     ARD_enviwrite_bands([v_input.l_dir,'/',n_map,'/CoverMap'],CoverMap,'uint8','bsq',all_yrs);
% %     clear CoverMap
% %     ARD_enviwrite_bands([v_input.l_dir,'/',n_map,'/CoverQAMap'],CoverQAMap,'uint8','bsq',all_yrs);
% %     clear CoverQAMap

end
% %% Add change event 
% change_map = double(enviread('CCDCMap/ChangeMap'));
% % minimum mapping unit
% num_obj = 5;
% 
% % spatial filtering
% for i = 1:size(change_map,3)
%     i
%     tmp = change_map(:,:,i);
%     tmp(tmp > 0) = 1;
%     segm_tmp=bwlabeln(tmp,8);
%     L = segm_tmp;
%     s = regionprops(L,'area');
%     area = [s.Area];
%     
%     % filter out cloud object < than num_cldoj pixels
%     idx = find(area >= num_obj);
%     tmp(ismember(L,idx)==0) = 0;
%     copy_map = change_map(:,:,i);
%     copy_map(tmp==0) = 0;
%     change_map(:,:,i) = copy_map;
% end
% % write filtered change map
% ARD_enviwrite_bands([v_input.l_dir,'/',n_map,'/FChangeMap'],change_map,'uint16','bsq',all_yrs,'example_img');
% 
% %% change analysis
% max_n = 30;
% map = double(enviread('CCDCMap/CoverMap1_4'));
% change_stat = [];
% 
% for i = 1:max_n-1
%     i
%     change_map = map(:,:,i)*100 + map(:,:,i+1);
%     idstb = map(:,:,i+1)-map(:,:,i);
%     change_stat = [change_stat;change_map(idstb~=0)];
% end
% 
% % update class number
% all_class = unique(change_stat);
% % update number of class
% n_class = length(all_class);
% % calculate proportion based # for each class
% number = hist(change_stat,all_class);
% prct = number/sum(number);

% %% add disturbance (3) to the classification map
% for i = 1:max_n
%     change_tmp = change_map(:,:,i);
%     class_tmp = cover_map(:,:,i);
%     class_tmp(change_tmp>0) = 3;
%     cover_map(:,:,i) = class_tmp;
% end



