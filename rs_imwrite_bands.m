function rs_imwrite_bands(varargin)
% rs_imwrite - Write an ENVI image with differnt bands with names (Copy Right - Zhe Zhu, Sept. 4th 2018)
% Examples
% rs_imwrite(image,'filename',info,n_names); 
% Input
%   image: image matrix
%   filename: File names
%   info: Image properties and cartographic information for ENVI or GeoTIFF file
%   n_names   : Band names 1985:2015 or 'disturb type';

if nargin < 4
    fprintf('Failed! Less Than Three Variables! \n');
    return
elseif nargin == 4    
    % Create input variables
    image = varargin{1};
    fname = varargin{2};
    info = varargin{3};
    n_names = varargin{4};
    
    im_size=size(image);
    im_size(3)=size(image,3);
    
    % ENVI data format numbers
    d = [1 2 3 4 5 12];
    
    % Check user input
    if ~ischar(fname)
        error('fname should be a char string');
    end
    
    cl=class(image);
    switch cl
        case 'uint8'
            t = d(1);
        case 'int16'
            t = d(2);
        case 'int32'
            t = d(3);
        case 'single'
            t = d(4);
        case 'double'
            t = d(5);
        case 'uint16'
            t = d(6);
        otherwise
            error('Data type not recognized');
    end
    
    % update info
    info.samples = im_size(2);
    info.lines = im_size(1);
    info.bands = im_size(3);
    info.data_type = t;
    
    try multibandwrite(image,fname,info.interleave);
        fprintf('Write ENVI Images %s',fname);
        fprintf((' . '));
    catch
        fprintf('Images Failed! Not ENVI Variables! \n');
        return
    end
    
    % Write header file
    fprintf((' . '));
    fid = fopen(strcat(fname,'.hdr'),'w');
    
    fprintf(fid,'%s \n','ENVI');
    n_field = fieldnames(info);
    n_value = struct2cell(info);
    
    for i = 1:length(n_field)-1
        n_field{i} = strrep(n_field{i},'_',' ');
        if ischar(n_value{i})
            fprintf(fid,'%s = %s \n',n_field{i},n_value{i});
        else
            fprintf(fid,'%s = %d \n',n_field{i},n_value{i});
        end
    end
    
    % band names
    fprintf(fid,'band names = {');
    if ischar(n_names)
        fprintf(fid,'%s',n_names);
    else
        for i=1:length(n_names)
            str_n = num2str(n_names(i));
            fprintf(fid,'%s',str_n);
            if i < length(n_names)
                fprintf(fid,',');
            end
        end
    end
    fprintf(fid,'}\n');
    
    %     elements={'samples =','lines   =','bands   =',...
%         'header offset =','file type =', 'data type =',...
%         'interleave =','sensor type =','byte order =',...
%         'map info =','coordinate system string =','wavelength units =',...
%         'band names ='};
    
%     fprintf(fid,'%s \n','ENVI');
%     fprintf(fid,'%s %s \n','description =',info.description);
%     fprintf(fid,'%s %d \n',elements{1,1},im_size(2));
%     fprintf(fid,'%s %d \n',elements{1,2},im_size(1));
%     fprintf(fid,'%s %d \n',elements{1,3},im_size(3));
%     fprintf(fid,'%s %d \n',elements{1,4},info.header_offset);
%     fprintf(fid,'%s %s \n',elements{1,5},info.file_type);
%     fprintf(fid,'%s %d \n',elements{1,6},t);
%     fprintf(fid,'%s %s \n',elements{1,7},info.interleave);
%     fprintf(fid,'%s %s \n',elements{1,8},info.sensor_type);
%     fprintf(fid,'%s %d \n',elements{1,9},info.byte_order);
%     fprintf(fid,'%s %s \n',elements{1,10},info.map_info);
%     fprintf(fid,'%s %s \n',elements{1,11},info.coordinate_system_string);
%     fprintf(fid,'%s %s \n',elements{1,12},info.wavelength_units);
    % fprintf(fid,'%s %s \n',elements{1,13},info.band_names);
    
    fclose(fid);
    fprintf((' . \n'));
else
    fprintf('Failed! Too Many Variables! \n');
    return
end