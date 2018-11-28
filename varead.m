function va=varead(filename,vaname)
% Read the value of the first vaname variable
% [im,jiDim,jiUL,resolu,ZC]=enviread(filename)
% Input file directory
% Output 1) im = data 
% Output 2) jiDim dimension cols and rows 
% Output 3)jiUL UpperLeft coner of the pixel (x,y)
% Output  4) resolution [x,y]
% Output 5) ZC zone code

% version 1.1 Able to read in USGS ENVI format data


fid_in1=fopen(filename);

if fid_in1~=-1
    fid_in=fid_in1;
else
    fprintf('Wrong file!\n');
    return;
end

geo_char=fscanf(fid_in,'%c',inf);
fclose(fid_in);

geo_char=geo_char';
geo_str=strread(geo_char,'%s');


va_id = strmatch(vaname,geo_str)+2;
va = char(geo_str(va_id(1)));
