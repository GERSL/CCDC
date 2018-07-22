function [jiDim,jiUL,resolu,ZC,bands]=envihdrread(filename)
% any ENVI file uint 8, int16, and uint16
% [jiDim,jiUL,resolu,ZC]=envihdrread(filename)
% Input file directory
% Output 1) im = data 
% Output 2) jiDim dimension cols and rows 
% Output 3)jiUL UpperLeft coner of the pixel (x,y)
% Output  4) resolution [x,y]
% Output 5) ZC zone code

filename_HDR=[filename,'.HDR'];
filename_hdr=[filename,'.hdr'];

fid_in1=fopen(filename_hdr,'r');
fid_in2=fopen(filename_HDR,'r');

if fid_in1~=-1
    fid_in=fid_in1;
elseif fid_in2~=-1
    fid_in=fid_in2;
else
    fprintf('Wrong ENVI header file!\n');
    fprintf('%s\n',filename); % show envi hdr file
    return;
end

geo_char=fscanf(fid_in,'%c',inf);
fclose(fid_in);

geo_char=geo_char';
geo_str=strread(geo_char,'%s');


indx_samples = strmatch('samples',geo_str)+2;
indx_lines = strmatch('lines',geo_str)+2;
indx_bands = strmatch('bands',geo_str)+2;
indx_datatype = strmatch('data',geo_str)+3;
indx_interleave = strmatch('interleave',geo_str)+2;
indx_xUL = strmatch('map',geo_str)+6;
indx_yUL = strmatch('map',geo_str)+7;
indx_xreso = strmatch('map',geo_str)+8;
indx_yreso = strmatch('map',geo_str)+9;
indx_zc = strmatch('map',geo_str)+10;

% read input image hdr
cols = str2double(geo_str(indx_samples));
rows = str2double(geo_str(indx_lines));
jiDim = [cols,rows];

bands = str2double(geo_str(indx_bands)); 
datatype = str2double(geo_str(indx_datatype));
interleave = char(geo_str(indx_interleave));
jiUL(1)=str2double(geo_str(indx_xUL)); 
jiUL(2)=str2double(geo_str(indx_yUL)); 
resolu(1)=str2double(geo_str(indx_xreso)); 
resolu(2)=str2double(geo_str(indx_yreso)); 
ZC=str2double(geo_str(indx_zc)); 

if datatype == 1
    in_type = 'uint8';
elseif datatype == 2
    in_type = 'int16';
elseif datatype == 3
    in_type = 'int32';
elseif datatype == 4
    in_type = 'single';
elseif datatype == 12
    in_type = 'uint16';
else
    error('Invalid read data type!');
end

% % Define the data set.
% % Read every other band of the data using the Band-Sequential format.
% im = multibandread(filename, [rows cols bands], ...
%                     [in_type,'=>',in_type], 0, interleave, 'ieee-le');
end
           
