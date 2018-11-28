function autoShowSynAll(varargin)
%AUTOSHOWSYNALL provde all synthetic data correpsonding to the prepared
%images.
% 
% Specific parameters
% ------------------------
%   'InDir'     Directory of CCDC change detection results.  
%                        Default is the path to current folder.
%   'OutDir'    Directory of the model output.  Default is the 
%                        same as InputDirectory.
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

parse(p,varargin{:});
dir_cur = p.Results.InDir;
dir_out = p.Results.OutDir;

% get num of total folders start with "L"
imf=dir(fullfile(dir_cur,'L*')); % folder names
% filter for Landsat folders
%imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
% WA timber harvest scene
imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
imf = [imf{:}];
imf = vertcat(imf{:});
% sort according to yeardoy
yeardoy = str2num(imf(:, 10:16));
[~, sort_order] = sort(yeardoy);
imf = imf(sort_order, :);

% number of images
num_t = size(imf,1);

for i = 1:num_t
    fprintf('Processing the %d th image.\r\n',i);
    % Find date for folder imf(i)
    yr = str2num(imf(i, 10:13));
    doy = str2num(imf(i, 14:16)); 
    autoShowSyn1('InputDirectory',dir_cur,...
        'OutputDirectory',dir_out,....
        'Year',yr,...
        'DOY',doy);
end
