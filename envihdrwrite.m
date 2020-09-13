function envihdrwrite(filename,dims,datatype,resolu,UL,interleave,zone)
%  write envi format file of unit8 and int16 
%  for example: enviwrite('testzz',data,'int16',[30,30],jiUL,'bsq',ZC,'this is for test only');
%  i,j => x, y => rows and cols
%  Input  1)  'filename' is the name of the file
%  Input  2)  'dims' is the dimenstion of the image (rows,cols,bands)
%  Input  3)  'dataype' is the datatype to write 'unit8' or 'int16'
%  Input  4)  'resolu' is the resolution of the pixel
%  Input  5)  'UL' is the UpperLeftPoint X Y of the UL pixel (not center of UL pixel)
%  Input  6)  'interleave' is the bsq bip bil format
%  Input  7)  'zone' is the UTM zone


% if strcmp(datatype,'uint8')
% %     envi_data=uint8(envi_data);
%     dt=1;
% elseif strcmp(datatype,'int16')
% %     envi_data=int16(envi_data);
%     dt=2;
% elseif strcmp(datatype,'uint16')
% %     envi_data=uint16(envi_data);
%     dt=12;
% else
%     fprintf('Invalid write data type!\n');
%     return;
% end

if strcmp(datatype,'uint8')
%    envi_data=uint8(envi_data);
    dt=1;
elseif strcmp(datatype,'int16')
%    envi_data=int16(envi_data);
    dt=2;
elseif strcmp(datatype,'int32')
%    envi_data=int32(envi_data);
    dt=3;
elseif strcmp(datatype,'uint16')
%    envi_data=uint16(envi_data);
    dt=12;
elseif strcmp(datatype,'single')
%    envi_data=single(envi_data);
    dt=4;
else
    fprintf('Invalid write data type!\n');
    return;
end

% n_dims=size(envi_data);

nrows=dims(1); 
ncols=dims(2); 
bands=1;

if length(dims)==3
    bands=dims(3);
end

% multibandwrite(envi_data,filename,interleave);

filename_hdr=[filename,'.hdr'];
                                
fid_out=fopen(filename_hdr,'wt');

fprintf(fid_out,'ENVI\n');
fprintf(fid_out,'description = {Landsat Scientific Data}\n');
fprintf(fid_out,'samples = %d\n',ncols); % samples is for j
fprintf(fid_out,'lines   = %d\n',nrows); % lines is for i
fprintf(fid_out,'bands   = %d\n',bands);
fprintf(fid_out,'header offset = 0\n');
fprintf(fid_out,'file type = ENVI Standard\n');
fprintf(fid_out,'data type = %d\n',dt);
fprintf(fid_out,'interleave = %s\n',interleave);
fprintf(fid_out,'sensor type = Landsat\n');
fprintf(fid_out,'byte order = 0\n');

if isempty(zone)
    %label as Landsat ARD
    fprintf(fid_out, 'map info = {Albers Conical Equal Area, 1, 1, %d, %d, %d, %d, WGS-84}\n',UL(1),UL(2),resolu(1),resolu(2));
    fprintf(fid_out, 'projection info = {9, 6378140, 6356755.288157528, 23, -96, 0, 0, 29.5, 45.5,WGS-84, Albers Conical Equal Area}\n');
    fprintf(fid_out, 'coordinate system string = {PROJCS["WGS_1984_Albers",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378140.0,298.257]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-96.0],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_origin",23.0],UNIT["Meter",1.0]]}');
else
    if (zone > 0)
        fprintf(fid_out, 'map info = {UTM, 1.000, 1.000, %d, %d, %d, %d, %d, North, WGS-84, units=Meters}',UL(1),UL(2),resolu(1),resolu(2),zone);
    else
        fprintf(fid_out, 'map info = {UTM, 1.000, 1.000, %d, %d, %d, %d, %d, South, WGS-84, units=Meters}',UL(1),UL(2),resolu(1),resolu(2),-zone);
    end
end

fclose(fid_out);
end



