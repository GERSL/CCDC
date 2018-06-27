function autoPrepareDataESPA(varargin)
%AUTOPREPAREDATAARD Prepare Landsat Surface Reflectance into CCDC format, 
% which are downloaded from USGS Earth Resources Observation and Science
% (EROS) Center Science Processing Architecture (ESPA)
% (https://espa.cr.usgs.gov/)
%
%   AUTOPREPAREDATAARD() automatically prepares all Landsat ESPA product in
%   the current folder into CCDC format.
%   AUTOPREPAREDATAARD(PARAM1,VAL1,PARAM2,VAL2,PARAM3,VAL3,PARAM4,VAL4) specifies 
%   parameters that control input and outout directory, the clear pixel
%   filter condition, and sample file (used to restrict same extent using 
%   nearest method based on the Topotoolbox 
%   https://topotoolbox.wordpress.com/topotoolbox/).
%
% Data Support
%   -------------
%   The input data must be of class .tif or .img. Only .tif format can be 
%   resamlped to a same extent.
%   
%
% Specific parameters
% ------------------------
%   'InputDirectory'     Directory of input data.  Default is the path to
%                        the current folder.
%   'OutputDirectory'    Directory of output data.  Default is the path to
%                        the current folder.
%   'ExtentSample'       An example geotiff file, of which extent will be 
%                        used as basic reference. All images will be 
%                        resampled to this same extent.  Default is not to
%                        do this process if no input for this.
%   'ClearPixelPercent'  Percentage of mininum clear pixels. Unit is %.
%                        Default is '20'.
%
%
%   Author:  Zhe Zhu (zhe.zhu#ttu.edu)
%            Shi Qiu (shi.qiu#ttu.edu)
%            Junxue Zhang (junxue.zhang#ttu.edu)
%   Date: 24. Jun, 2018

    %% get parameters from inputs
    % where the all Landsat zipped files are
    dir_cur = pwd;
    % where the output files are
    dir_out = dir_cur;
    % min clear pixel
    clr_pct_min = 20; % unit %
    % total number of bands
    nbands = 8;
    
    p = inputParser;
    p.FunctionName = 'prepParas';
    % optional
    % default values.
    addParameter(p,'InputDirectory',dir_cur);
    addParameter(p,'OutputDirectory',dir_out);
    addParameter(p,'ClearPixelPercent',clr_pct_min);
    addParameter(p,'ExtentSample','');
    
    % request user's input
    parse(p,varargin{:});
    dir_cur = p.Results.InputDirectory;
    dir_out = p.Results.OutputDirectory;
    clr_pct_min = p.Results.ClearPixelPercent;
    trgt_file = p.Results.ExtentSample;
    

    %% Locate to the current directory
    % name of the temporary folder for extracting zip files
    name_tmp = 'tmp';
    % remove all temp folders
    tmpf = dir(fullfile(dir_out,[name_tmp,'*']));
    if ~isempty(tmpf)
        for i = 1:length(tmpf)
            if tmpf(i).isdir
                rmdir(fullfile(dir_out,tmpf(i).name),'s');
            end
        end
    end
    
    %% Filter for Landsat folders
    % get num of total folders start with "L"
    %imf = dir('L*SR.tar'); % folder names
    imfs = dir(fullfile(dir_cur,'L*.tar.gz'));
    % filter for Landsat folders
    % espa data
    imfs = regexpi({imfs.name}, 'L(T05|T04|E07|C08)(\w*)\-(\w*).tar.gz', 'match'); 
    imfs = [imfs{:}];
    if isempty(imfs)
        warning('No images here!');
        return;
    end
    imfs = vertcat(imfs{:});
    % sort according to yeardoy
    yyyymmdd = str2num(imfs(:, 11:18)); % should change for different sets
    [~, sort_order] = sort(yyyymmdd);
    imfs = imfs(sort_order, :);
    % number of folders start with "L"
    num_t = size(imfs,1);
    fprintf('A total of %d images will be prepared...\n',num_t);
    for i = 1:num_t
        % name of the temporary folder for extracting zip files
        n_tmp = [name_tmp,num2str(i)];
        imf = imfs(i,:);
        % new filename in format of LXSPPPRRRYYYYDOYLLLTT
        yr = str2num(imf(11:14));
        % converst mmdd to doy
        mm = str2num(imf(15:16));
        dd = str2num(imf(17:18));
        doy = datenummx(yr,mm,dd)-datenummx(yr,1,0);
        % set folder and image name
        n_mtl = [imf([1,2,4:14]),num2str(doy,'%03d'),'0',imf(19:22)];
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
                fprintf('%s exsit in stacked images folder\n',tmp_img);
                continue;
            end
        end
        % unzip the images to a temporary folder
        %fprintf('Unzip the %dth image ...\n',i);
        try
            n_gun = gunzip(fullfile(dir_cur,imf),fullfile(dir_out,n_tmp));
            n_tar = untar(fullfile(dir_out,n_tmp,imf(1:end-3)),fullfile(dir_out,n_tmp));
        catch me
            if isfolder(fullfile(dir_out,n_tmp))
                rmdir(fullfile(dir_out,n_tmp),'s');
            end
            fprintf('File cannot be found in the %dth image',i);
            continue;
        end

        % decide image format (tif or envi)
        env_cfmask = dir(fullfile(dir_out,n_tmp,'L*pixel_qa.img'));
        tif_cfmask = dir(fullfile(dir_out,n_tmp,'L*pixel_qa.tif'));

        % picking surf ref 1-7, bt, and cfmask and save to the images folder
        % get names of surf 1-7
        % fprintf('Reading images ...\n');

        % read cfmask first to caculate clear pixel percet
        if ~isempty(env_cfmask) % envi format
            env_cfmask = fullfile(dir_out,n_tmp,env_cfmask.name);
            [cfmask0,jidim,jiul,resolu,zc] = enviread(env_cfmask);
        else
            tif_cfmask = fullfile(dir_out,n_tmp,tif_cfmask.name);
            cfmask0 = geotiffread(tif_cfmask);
        end
        cfmask = cfmask0;

        % convert pixel QA to fmask values
        cfmask(bitget(cfmask0,1) == 1) = 255;
        cfmask(bitget(cfmask0,2) == 1) = 0;
        cfmask(bitget(cfmask0,3) == 1) = 1;
        cfmask(bitget(cfmask0,4) == 1) = 2;
        cfmask(bitget(cfmask0,5) == 1) = 3;
        cfmask(bitget(cfmask0,6) == 1) = 4;

        clr_pct = sum(cfmask(:)<=1)/sum(cfmask(:)<255);
        clr_pct = 100*clr_pct;
        if clr_pct < clr_pct_min % less than 20% clear observations
            % remove the tmp folder
            % fprintf('Clear observation less than 20 percent (%.2f) ...\n',clr_pct*100);
            rmdir(fullfile(dir_out,n_tmp),'s');
            fprintf('Clear observation less than %.2f percent (%.2f) ...\n',clr_pct_min,clr_pct*100);
            continue;
        else
            if ~isempty(tif_cfmask) % tif format(tif_cfmask);
                % get projection information from geotiffinfo
                info = geotiffinfo(tif_cfmask);
                jidim = [info.SpatialRef.RasterSize(2),info.SpatialRef.RasterSize(1)];
                jiul = [info.SpatialRef.XLimWorld(1),info.SpatialRef.YLimWorld(2)];
                resolu = [info.PixelScale(1),info.PixelScale(2)];
                zc = info.Zone;
            end

            % prelocate image for the stacked image
            stack = zeros(jidim(2),jidim(1),nbands,'int16');
            % give cfmask to the last band
            stack(:,:,end) = cfmask;
        end

        if ~isempty(env_cfmask) % envi format
            if str2num(n_mtl(3)) < 8
                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band1.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b1 = enviread(n_surf);
                stack(:,:,1) = surf_b1;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band2.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b2 = enviread(n_surf);
                stack(:,:,2) = surf_b2;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band3.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b3 = enviread(n_surf);
                stack(:,:,3) = surf_b3;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band4.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b4 = enviread(n_surf);
                stack(:,:,4) = surf_b4;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band5.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b5 = enviread(n_surf);
                stack(:,:,5) = surf_b5;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band7.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b7 = enviread(n_surf);
                stack(:,:,6) = surf_b7;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*bt_band6.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b6 = enviread(n_surf);
                stack(:,:,7) = surf_b6;

            else

                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band2.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b1 = enviread(n_surf);
                stack(:,:,1) = surf_b1;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band3.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b2 = enviread(n_surf);
                stack(:,:,2) = surf_b2;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band4.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b3 = enviread(n_surf);
                stack(:,:,3) = surf_b3;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band5.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b4 = enviread(n_surf);
                stack(:,:,4) = surf_b4;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band6.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b5 = enviread(n_surf);
                stack(:,:,5) = surf_b5;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band7.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b7 = enviread(n_surf);
                stack(:,:,6) = surf_b7;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*bt_band10.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b6 = enviread(n_surf);
                stack(:,:,7) = surf_b6;
            end
        else
            if isempty(trgt_file) % no sample file as target
                if str2num(n_mtl(3)) < 8
                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band1.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b1 = geotiffread(n_surf);
                    stack(:,:,1) = surf_b1;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band2.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b2 = geotiffread(n_surf);
                    stack(:,:,2) = surf_b2;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band3.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b3 = geotiffread(n_surf);
                    stack(:,:,3) = surf_b3;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band4.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b4 = geotiffread(n_surf);
                    stack(:,:,4) = surf_b4;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band5.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b5 = geotiffread(n_surf);
                    stack(:,:,5) = surf_b5;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band7.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b7 = geotiffread(n_surf);
                    stack(:,:,6) = surf_b7;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*bt_band6.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b6 = geotiffread(n_surf);
                    stack(:,:,7) = surf_b6;

                else

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band2.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b1 = geotiffread(n_surf);
                    stack(:,:,1) = surf_b1;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band3.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b2 = geotiffread(n_surf);
                    stack(:,:,2) = surf_b2;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band4.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b3 = geotiffread(n_surf);
                    stack(:,:,3) = surf_b3;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band5.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b4 = geotiffread(n_surf);
                    stack(:,:,4) = surf_b4;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band6.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b5 = geotiffread(n_surf);
                    stack(:,:,5) = surf_b5;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band7.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b7 = geotiffread(n_surf);
                    stack(:,:,6) = surf_b7;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*bt_band10.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b6 = geotiffread(n_surf);
                    stack(:,:,7) = surf_b6;
                end
            else
                % have targt file
                try
                    trgt_obj = GRIDobj(trgt_file);
                    trgt_obj.Z = [];% empty memory
                catch
                    fprintf('Sample file can not be support. Only geotiff is workable.\n');
                    return;
                end
                
                if str2num(n_mtl(3)) < 8
                    
                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band1.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b1 = geotiffread(n_surf);
                    % create GRIDobj
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b1; clear surf_b1;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,1) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band2.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b2 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b2; clear surf_b2;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,2) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band3.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b3 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b3; clear surf_b3;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,3) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band4.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b4 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b4; clear surf_b4;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,4) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band5.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b5 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b5; clear surf_b5;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,5) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band7.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b7 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b7; clear surf_b7;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,6) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*bt_band6.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b6 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b6; clear surf_b6;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,7) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                else

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band2.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b1 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b1; clear surf_b1;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,1) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band3.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b2 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b2; clear surf_b2;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,2) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band4.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b3 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b3; clear surf_b3;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,3) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band5.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b4 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b4; clear surf_b4;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,4) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band6.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b5 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b5; clear surf_b5;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,5) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*sr_band7.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b7 = geotiffread(n_surf);=
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b7; clear surf_b7;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,6) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;

                    n_surf = dir(fullfile(dir_out,n_tmp,'L*bt_band10.tif'));
                    n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                    surf_b6 = geotiffread(n_surf);
                    tmp_obj = GRIDobj(n_surf);
                    % gridobj reader may be NaN because this is special for DEM data.
                    tmp_obj.Z = surf_b6; clear surf_b61;
                    % same extent and resolution
                    tmp_obj_same_extn = resample(tmp_obj,trgt_obj,'nearest',true);
                    stack(:,:,7) = tmp_obj_same_extn.Z;
                    clear tmp_obj tmp_obj_same_extn;
                end
            end
        end

        % new stacked bip image
        % get image name from cfmask file
        %     n_mtl = dir([dir_cur,n_tmp,'L*_pixel_qa.tif']);
        %     n_mtl = strsplit(char(n_mtl.name),'_'); % remove .txt
        %     n_mtl = strjoin(n_mtl([1 3 4 6 7]),''); % [1 3 4 7 8]for Espa files

        n_stack = [char(n_mtl),'_MTLstack'];

        % generate folder name
        %        n_mtl = strsplit(char(n_mtl),'_'); % remove _MTL
        %        n_mtl = char(n_mtl(1));

        % add directory
        n_dir = fullfile(dir_out,n_mtl);
        if fold_exist == 0
            mkdir(n_dir);
        end

        % write to images folder
        % fprintf('Writing %s image ...\n',n_mtl);
        n_stack = fullfile(n_dir,n_stack);
        enviwrite(n_stack,stack,'int16',resolu,jiul,'bip',zc);

        % remove the tmp folder
        rmdir(fullfile(dir_out,n_tmp),'s');
    end
end
