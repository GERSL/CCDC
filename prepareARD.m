function [outputArg1,outputArg2] = prepareARD(dir_cur,dir_out,imf)
%PREPAREARD prepare single Landsat ARD.
% Inputs:
% dir_cur: directory of input
% dir_out: directory of output
% imf: name of image without any suffix


% Outputs:
% info = 0: 
% info = 1: BT or SR file cannot be found
        info = 0;
        % new filename in format of LXSPPPRRRYYYYDOYLLLTT
        yr = str2num(imf(16:19));
        % converst mmdd to doy
        mm = str2num(imf(20:21));
        dd = str2num(imf(22:23));
        doy = datenummx(yr,mm,dd)-datenummx(yr,1,0);
        % set folder and image name
        n_mtl = [imf([1,2,4,9:14,16:19]),num2str(doy,'%03d'),imf([34,16,38:40])];
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
            sr_tar = untar(fullfile(dir_cur,imf),fullfile(dir_out,n_tmp));
            if ispc
                % Code to run on Windows platform
                bt_tar = untar(fullfile(dir_cur,[imf(1:end-6),'BT.tar']),fullfile(dir_out,n_tmp));
            else
                bt_tar = untar(fullfile(dir_cur,[imf(1:end-2),'BT']),fullfile(dir_out,n_tmp));
            end
        catch me
            if isfolder(fullfile(dir_out,n_tmp))
                rmdir(fullfile(dir_out,n_tmp),'s');
            end
            fprintf('BT or SR file cannot be found in the %dth image',i);
            continue;
        end

        % decide image format (tif or envi)
        env_cfmask = dir(fullfile(dir_out,n_tmp,'L*PIXELQA.img'));
        tif_cfmask = dir(fullfile(dir_out,n_tmp,'L*PIXELQA.tif'));

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
            % fprintf('Clear pixels less than %.2f percent (%.2f) ...\n',clr_pct_min,clr_pct);
            fprintf('Clear pixels less than %.2f percent (%.2f) for %s\n',clr_pct_min,clr_pct,imf);
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

        if ~isempty(env_cfmask) % envi forma
            if str2num(n_mtl(3)) < 8
                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB1.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b1 = enviread(n_surf);
                stack(:,:,1) = surf_b1;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB2.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b2 = enviread(n_surf);
                stack(:,:,2) = surf_b2;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB3.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b3 = enviread(n_surf);
                stack(:,:,3) = surf_b3;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB4.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b4 = enviread(n_surf);
                stack(:,:,4) = surf_b4;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB5.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b5 = enviread(n_surf);
                stack(:,:,5) = surf_b5;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB7.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b7 = enviread(n_surf);
                stack(:,:,6) = surf_b7;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*BTB6.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b6 = enviread(n_surf);
                stack(:,:,7) = surf_b6;

            else

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB2.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b1 = enviread(n_surf);
                stack(:,:,1) = surf_b1;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB3.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b2 = enviread(n_surf);
                stack(:,:,2) = surf_b2;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB4.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b3 = enviread(n_surf);
                stack(:,:,3) = surf_b3;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB5.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b4 = enviread(n_surf);
                stack(:,:,4) = surf_b4;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB6.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b5 = enviread(n_surf);
                stack(:,:,5) = surf_b5;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB7.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b7 = enviread(n_surf);
                stack(:,:,6) = surf_b7;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*BTB10.img'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b6 = enviread(n_surf);
                stack(:,:,7) = surf_b6;
            end
        else
            if str2num(n_mtl(3)) < 8

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB1.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b1 = geotiffread(n_surf);
                stack(:,:,1) = surf_b1;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB2.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b2 = geotiffread(n_surf);
                stack(:,:,2) = surf_b2;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB3.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b3 = geotiffread(n_surf);
                stack(:,:,3) = surf_b3;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB4.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b4 = geotiffread(n_surf);
                stack(:,:,4) = surf_b4;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB5.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b5 = geotiffread(n_surf);
                stack(:,:,5) = surf_b5;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB7.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b7 = geotiffread(n_surf);
                stack(:,:,6) = surf_b7;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*BTB6.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b6 = geotiffread(n_surf);
                stack(:,:,7) = surf_b6;

            else

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB2.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b1 = geotiffread(n_surf);
                stack(:,:,1) = surf_b1;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB3.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b2 = geotiffread(n_surf);
                stack(:,:,2) = surf_b2;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB4.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b3 = geotiffread(n_surf);
                stack(:,:,3) = surf_b3;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB5.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b4 = geotiffread(n_surf);
                stack(:,:,4) = surf_b4;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB6.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b5 = geotiffread(n_surf);
                stack(:,:,5) = surf_b5;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*SRB7.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b7 = geotiffread(n_surf);
                stack(:,:,6) = surf_b7;

                n_surf = dir(fullfile(dir_out,n_tmp,'L*BTB10.tif'));
                n_surf = fullfile(dir_out,n_tmp,n_surf.name);
                surf_b6 = geotiffread(n_surf);
                stack(:,:,7) = surf_b6;
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
        % if ~isfolder(n_dir)
        %     mkdir(n_dir);
        % end
        % works before 2017b
        n_dir_isfolder = dir(n_dir);
        if isempty(n_dir_isfolder)
            mkdir(n_dir);
        end
        clear n_dir_isfolder;

        % write to images folder
        % fprintf('Writing %s image ...\n',n_mtl);
        n_stack = fullfile(n_dir,n_stack);
        enviwrite(n_stack,stack,'int16',resolu,jiul,'bip',zc);
end

