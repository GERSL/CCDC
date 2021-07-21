function autoPrepareDataESPAC2(varargin)
%AUTOPREPAREDATAARDC2 Prepare Landsat Collection 2 Surface Reflectance into CCDC format, 
% which are downloaded from USGS Earth Resources Observation and Science
% (EROS) Center Science Processing Architecture (ESPA)
% (https://espa.cr.usgs.gov/)
%
%   AUTOPREPAREDATAARDC2() automatically prepares all Landsat ESPA product in
%   the current folder into CCDC format.
%   AUTOPREPAREDATAARDC2(PARAM1,VAL1,PARAM2,VAL2,PARAM3,VAL3,PARAM4,VAL4) specifies 
%   parameters that control input and outout directory, the clear pixel
%   filter condition, and sample file (used to restrict same extent using 
%   nearest method based on the Topotoolbox 
%   https://topotoolbox.wordpress.com/topotoolbox/).
%
% Data Support
%   -------------
%   The input data must be of class .tif.  Only .tif format can be 
%   resamlped to a same extent.
%   
%
% Specific parameters
% ------------------------
%   'pathin'          Directory of input Landsat Collection 2 Level 2 data.  Default is the path to
%                         the current folder.
%   'pathout'        Directory of outputing the stacked data.  Default is the path same as 'pathin'
%   'extfile'          An example geotiff file, of which extent will be used as basic reference. All images will be 
%                        resampled to this same extent, with nearest. If no specific, will stack the original Landsat data.
%   'clear'           Percentage of mininum clear pixels (non-ice/snow covered). Unit is %. Default is '20'.
%
%   task (optional):        Task ID of parallel computation
%
%   ntasks (optional):      Total number of tasks of parallel computation
%
%   'msg'            [true or false] Display process status message. Default is ture.
%
%
% Example:
%   autoPrepareDataESPAC2('pathin' , directory, 'extfile' , filepath, 'clear', 0);
%
%   Author:  
%               Shi Qiu (shi.qiu#uconn.edu) & Zhe Zhu (zhe#uconn.edu)
%               GERS Lab, UCONN
%                 
%
%   This version was tested on Matlab 2021b
%   Date: 20 July, 2021

    warning('off','all'); % do not show warning information
    % not to display the conflict between GRIDOBJ and internal funcitons
    % not to display the different projection
    
    %% Add search path
    addpath(fullfile(fileparts(mfilename('fullpath')), 'GRIDobj'));

    %% Parameters from inputs
    p = inputParser;
    p.FunctionName = 'prepParas';
    % optional
    % default values.
    addParameter(p,'pathin', pwd);
    addParameter(p,'pathout', '');
    addParameter(p,'clear', 20 ); % unit %
    addParameter(p,'extfile','');
    addParameter(p,'task', 1); % 1st task
    addParameter(p,'ntasks', 1); % single task to compute
    addParameter(p,'msg', true);
    
    % request user's input
    parse(p,varargin{:});
    dir_cur = p.Results.pathin;
    dir_out = p.Results.pathout;
    if isempty(dir_out)
        dir_out = dir_cur;
    end
    clr_pct_min = p.Results.clear;
    task = p.Results.task;
    ntasks = p.Results.ntasks;
    msg = p.Results.msg;
    trgt_file = p.Results.extfile;
    
    %% Locate to the current directory
    % name of the temporary folder for extracting zip files
    name_tmp = 'tmp_';
    
    
    %% Filter for Landsat folders
    % get num of total folders start with "L", with end of "SR.tar"
    imfs = dir(fullfile(dir_cur,'L*.tar'));
    % filter for Landsat folders
    % espa data
    imfs = regexpi({imfs.name}, 'L(T05|T04|E07|C08)(\w*)\_(\w*).tar', 'match'); 
    if isempty(imfs)
        warning( fprintf('No Landsat images can be found at %d! \n', dir_cur));
        return;
    end
    imfs = [imfs{:}];
    imfs = vertcat(imfs{:});
    % sort according to yeardoy
    yyyymmdd = str2num(imfs(:, 18:25)); % should change for different sets
    [~, sort_order] = sort(yyyymmdd);
    imfs = imfs(sort_order, :);
    clear yyyymmdd sort_order;
   
    
    % total number of bands
    nbands = 8; % i.e., Blue, Green, Red, NIR, SWIR1, SWIR2, TIR, and Fmask QA.
    
    
     %% Assign stacking tasks to each core
    num_t = size(imfs,1);
    tasks_per = ceil(num_t/ntasks);
    start_i = (task-1)*tasks_per + 1;
    end_i = min(task*tasks_per, num_t);
    if msg
        fprintf('At task# %d/%d to stack %d (of %d images)\n', task, ntasks, end_i- start_i + 1, num_t);
    end
    % create a task folder, that will be uesed to store the divided row
    % data for all the images @ the currenr core
    dir_out_tmp = fullfile(dir_out, sprintf('tmptaskfolder_%d_%d', task, ntasks));
    if ~isfolder(dir_out_tmp)
        mkdir(dir_out_tmp);
    end
    
    
    for i = start_i:end_i
        % name of Landsat image
        imf = imfs(i,:);
        % name of the temporary folder for extracting zip files
        n_tmp = [name_tmp, imf];
        % new filename in format of LXSPPPRRRYYYYDOYLLLTT, like LE70210372020181002T1
        % converst year, mm, dd, and doy
        yr = str2double(imf(18:21));
        mm = str2double(imf(22:23));
        dd = str2double(imf(24:25));
        doy = datenummx(yr, mm, dd)-datenummx(yr, 1, 0);
        % set folder and image name such as LE70210372020181002T1
        n_mtl = [imf([1, 2, 4,11:16]), ... % Sensor and Path and Row
            num2str(yr,'%04d'), ... % Year
            num2str(doy,'%03d'), ... % DOY
            '0', imf(19:20), ... % Collection 2
            imf(39:40)]; % T1 or T2
        % check if folder exsit or not
        % names of image folder that are processed
        n_img = dir(fullfile(dir_out,'L*'));
        num_img = size(n_img,1);
        % check all folders we have
        % record exist or not
        rec_exist = 0;
        if num_img > 0
            for i_check = 1:num_img
                if n_img(i_check).isdir
                    % each image folder name
                    tmp_img = n_img(i_check).name;
        
                    tmp_zip = n_mtl;
                    if strcmp(tmp_img(1:16),tmp_zip(1:16))
                        outf = dir(fullfile(dir_out,tmp_img,[char(n_mtl),'_MTLstack']));
                        if ~isempty(outf)
                            rec_exist = 1;
                            break;
                        end
                    end
                end
            end
            % continue if the folder already exist
            if rec_exist > 0
                    if msg
                        fprintf('%s exsit in stacked images folder\n',tmp_img);
                    end
                continue;
            end
        end
        % unzip the images to a temporary folder
        if msg
            fprintf('Untar the %dth image ...\n',i);
        end
        try
            untar(fullfile(dir_cur,imf),fullfile(dir_out_tmp,n_tmp));
        catch
            if isfolder(fullfile(dir_out_tmp,n_tmp))
                rmdir(fullfile(dir_out_tmp,n_tmp),'s');
            end
            fprintf('Error occurs in untaring the %dth image \n',i);
            continue;
        end
        imf = imf(: ,1:40); % remove .tar
        %% Cloud cover
        % read cfmask first to caculate clear pixel percet
        tif_cfmask = fullfile(dir_out_tmp, n_tmp, [imf, '_QA_PIXEL.tif']);
        cfmask0 = GRIDobj(tif_cfmask);
        cfmask = cfmask0.Z;
        % convert pixel QA to fmask values
        % see more details from USGS document at https://prd-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/atoms/files/LSDS-1328_Landsat8-9-OLI-TIRS-C2-L2-DFCB-v6.pdf
        cfmask(bitget(cfmask0.Z,1) == 1) = 255; % Filled
        cfmask(bitget(cfmask0.Z,7) == 1) = 0; % Clear Land and Water [No Cloud & No Dialted Cloud]
        cfmask(bitget(cfmask0.Z,8) == 1) = 1; % Water
        cfmask(bitget(cfmask0.Z,5) == 1) = 2; % Cloud Shadow
        cfmask(bitget(cfmask0.Z,6) == 1) = 3; % Snow
        cfmask(bitget(cfmask0.Z,2) == 1) = 4; % Dilated Cloud
        cfmask(bitget(cfmask0.Z,4) == 1) = 4; % Cloud

        clr_pct = sum(cfmask(:)<=1)/sum(cfmask(:)<255);
        clr_pct = 100*clr_pct;
        if clr_pct < clr_pct_min % less than 20% clear observations
            % remove the tmp folder
            % fprintf('Clear observation less than 20 percent (%.2f) ...\n',clr_pct*100);
            rmdir(fullfile(dir_out_tmp,n_tmp),'s');
            % fprintf('Clear pixels less than %.2f percent (%.2f) ...\n',clr_pct_min,clr_pct););
            if msg
                fprintf('Clear pixels less than %.2f percent (%.2f) for %s\n',clr_pct_min,clr_pct,imf);
            end
            continue;
        else
            %% Geo info
            if ~isempty(trgt_file)
                trgt_obj = GRIDobj(trgt_file);
                trgt_obj.Z = zeros(trgt_obj.size)+255;% empty memory initially masked as 255

                % update cfmask
                cfmask0.Z = cfmask;
                clear cfmask;
                % same extent and resolution
                tmp_obj_same_extn = resample(cfmask0,trgt_obj,'nearest',true,'fillval',255);
                cfmask = tmp_obj_same_extn.Z;
                clear tmp_obj_same_extn;
        
                % get projection information from geotiffinfo
                try % normal geotiff
                        info = geotiffinfo(trgt_file);
                        jidim = [info.SpatialRef.RasterSize(2),info.SpatialRef.RasterSize(1)];                    
                        try
                            jiul = [info.SpatialRef.XLimWorld(1),info.SpatialRef.YLimWorld(2)];
                        catch
                            jiul = [info.SpatialRef.LongitudeLimits(1),info.SpatialRef.LatitudeLimits(2)];
                        end
                        resolu = [info.PixelScale(1),info.PixelScale(2)];
                        zc = info.Zone;
                catch % cloud-format geotiff
                        % get projection information from geotiffinfo
                        info = georasterinfo(trgt_file);
                        jidim = [info.RasterSize(2), info.RasterSize(1)];                    
                        jiul = [info.RasterReference.XLimWorld(1), info.RasterReference.YLimWorld(2)];
                        resolu = [info.RasterReference.SampleSpacingInWorldX, info.RasterReference.SampleSpacingInWorldY];
                        % to have UTM Zone
                        zc = info.CoordinateReferenceSystem.Name;
                        zc = strsplit(zc, ' ');
                        zc = char(zc(end)); % i.e., 16N
                        if strcmpi( zc(end), 'n')
                            zc = str2double(zc(1:end-1));
                        else
                            zc = 0 - str2double(zc(1:end-1));
                        end
                end
                clear info;
            else
                % get projection information from geotiffinfo
                info = georasterinfo(tif_cfmask);
                jidim = [info.RasterSize(2), info.RasterSize(1)];                    
                jiul = [info.RasterReference.XLimWorld(1), info.RasterReference.YLimWorld(2)];
                resolu = [info.RasterReference.SampleSpacingInWorldX, info.RasterReference.SampleSpacingInWorldY];
                % to have UTM Zone
                zc = info.CoordinateReferenceSystem.Name;
                zc = strsplit(zc, ' ');
                zc = char(zc(end)); % i.e., 16N
                if strcmpi( zc(end), 'n')
                    zc = str2double(zc(1:end-1));
                else
                    zc = 0 - str2double(zc(1:end-1));
                end
                clear info;
            end
            clear cfmask0;
            %% stack each band
            % prelocate image for the stacked image
            stack = zeros(jidim(2),jidim(1),nbands,'int16');
            % give cfmask to the last band
            stack(:,:,end) = cfmask;
            if isempty(trgt_file) % no sample file as target
                if str2double(n_mtl(3)) < 8
                    stack(:,:,1) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B1.TIF']));
                    stack(:,:,2) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B2.TIF']));
                    stack(:,:,3) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B3.TIF']));
                    stack(:,:,4) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B4.TIF']));
                    stack(:,:,5) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B5.TIF']));
                    stack(:,:,6) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B7.TIF']));
                    stack(:,:,7) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_ST_B6.TIF']));
                else
                    stack(:,:,1) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B2.TIF']));
                    stack(:,:,2) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B3.TIF']));
                    stack(:,:,3) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B4.TIF']));
                    stack(:,:,4) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B5.TIF']));
                    stack(:,:,5) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B6.TIF']));
                    stack(:,:,6) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B7.TIF']));
                    stack(:,:,7) = readgeoraster(fullfile(dir_out_tmp,n_tmp, [imf, '_ST_B10.TIF']));
                end
            else
                if str2num(n_mtl(3)) < 8
                    stack(:,:,1) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B1.TIF']), trgt_obj);
                    stack(:,:,2) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B2.TIF']), trgt_obj);
                    stack(:,:,3) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B3.TIF']), trgt_obj);
                    stack(:,:,4) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B4.TIF']), trgt_obj);
                    stack(:,:,5) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B5.TIF']), trgt_obj);
                    stack(:,:,6) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B7.TIF']), trgt_obj);
                    stack(:,:,7) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_ST_B6.TIF']), trgt_obj);
                else
                    stack(:,:,1) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B2.TIF']), trgt_obj);
                    stack(:,:,2) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B3.TIF']), trgt_obj);
                    stack(:,:,3) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B4.TIF']), trgt_obj);
                    stack(:,:,4) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B5.TIF']), trgt_obj);
                    stack(:,:,5) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B6.TIF']), trgt_obj);
                    stack(:,:,6) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_SR_B7.TIF']), trgt_obj);
                    stack(:,:,7) = wrapimage(fullfile(dir_out_tmp,n_tmp, [imf, '_ST_B10.TIF']), trgt_obj);
                end
            end
        end

        % convert to the same format as we used before, and then we do
        % not need to change code of CCD.
        % examine the scale factors from .xml or 
        % https://prd-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/atoms/files/LSDS-1619_Landsat8-C2-L2-ScienceProductGuide-v2.pdf
        stack(:,:,1) = 10000.*(double(stack(:,:,1)).*0.0000275 - 0.2); % rescale to 0 ~10000
        stack(:,:,2) = 10000.*(double(stack(:,:,2)).*0.0000275 - 0.2); % rescale to 0 ~10000
        stack(:,:,3) = 10000.*(double(stack(:,:,3)).*0.0000275 - 0.2); % rescale to 0 ~10000
        stack(:,:,4) = 10000.*(double(stack(:,:,4)).*0.0000275 - 0.2); % rescale to 0 ~10000
        stack(:,:,5) = 10000.*(double(stack(:,:,5)).*0.0000275 - 0.2); % rescale to 0 ~10000
        stack(:,:,6) = 10000.*(double(stack(:,:,6)).*0.0000275 - 0.2); % rescale to 0 ~10000
        stack(:,:,7) =    10.*(double(stack(:,:,7)).*0.00341802 + 149); % rescale same as Landsat Collection 1 ARD (0.1)


        n_stack = [char(n_mtl),'_MTLstack'];
        
        % add directory
        n_dir = fullfile(dir_out,n_mtl);
        if ~isfolder(n_dir)
            mkdir(n_dir);
        end

        % write to images folder
%         fprintf('Writing %s image ...\n',n_mtl);
        n_stack = fullfile(n_dir,n_stack);
        enviwrite(n_stack,stack,'int16',resolu,jiul,'bip',zc);

        % remove the tmp folder
        rmdir(fullfile(dir_out_tmp,n_tmp),'s');
    end
    % remove the tmp parent folder
    rmdir(fullfile(dir_out_tmp),'s');
end

function stack = wrapimage(n_surf, trgt_obj)
        tmp_obj = GRIDobj(n_surf);
        % same extent and resolution
        tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true,'fillval',-9999);
        stack = tmp_obj_same_extn.Z;
        clear tmp_obj tmp_obj_same_extn;
end
