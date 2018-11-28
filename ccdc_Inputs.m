function v_input = ccdc_Inputs
% function v_input=Variable_Inputs
% Function for input variables and paths
% Use the default fmask toobox developed by Zhe Zhu
% addpath('~/ccdc');
% Tools of RFC
% addpath('~/ccdc/RF_Class_C'); % not good 
% addpath('~/Algorithms/CCDC/Tools/RF_Class_C');
%% Inputs:
% 1. Location of "im_roi" & place to save "modelRF"

% Canada beetle infestation project
% v_input.l_dir='/projectnb/landsat/projects/Beetle_Infest_Map/images/';

% LCMS project
% v_input.l_dir='/projectnb/landsat/projects/LCMS/4530/images/'; % done 
% v_input.l_dir='/projectnb/landsat/projects/LCMS/2727/images/'; % done 13
% v_input.l_dir='/projectnb/landsat/projects/LCMS/1228/images/'; % done
% v_input.l_dir='/projectnb/landsat/projects/LCMS/1432/images/'; % done
% v_input.l_dir='/projectnb/landsat/projects/LCMS/1637/images/'; % done
% v_input.l_dir='/projectnb/landsat/projects/LCMS/3532/images/'; % done 9

% Guangzhou project
% v_input.l_dir='/projectnb/landsat/users/fuyc/p122r044/images/'; 

% ACRE project
% v_input.l_dir='/projectnb/landsat/projects/ACRE/stacks/p004r066/images/';
% v_input.l_dir='/projectnb/landsat/projects/ACRE/stacks/p005r066/images/';
% v_input.l_dir='/projectnb/landsat/projects/ACRE/stacks/p003r067/images/';
% v_input.l_dir='/projectnb/landsat/projects/ACRE/stacks/p002r067/images/';



% LCMAP
% Five sites for USGS training data strategy development
% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p046r027/images'; %
% (1) 1.NLCD ANC done! 2.Fmask_stat done! 3.CCD re-done! 4.ChangeNum done!
%     5.Trends done! 6. Preparation re-done! 

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p043r034/images'; %
% (2) 1.NLCD ANC done! 2.Fmask_stat done! 3.CCD re-done! 4. ChangeNum done!
%     5.Trends done! 6. Preparation re-done!

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p028r033/images'; % 
% (3) 1.NLCD ANC done! 2.Fmask_stat done! 3.CCD re-done! 4.ChangeNum done!
%     5.Trends done! 6. Preparation re-done! 

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p023r037/images'; % 
% (4) 1.NLCD ANC done! 2.Fmask stat done! 3. CCD re-done! 4. ChangeNum done!
%     5.Trends done! 6. Preparation re-done!

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p027r027/images'; % 
% (5) 1. NLCD ANC done! 2.Fmask_stat done! 3.CCD re-done! 4.ChangeNum done!
%     5.Trends done! 6. Preparation re-done!

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p033r029/images'; % (6) 

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p016r040/images'; % (7) 

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p035r032/images'; % (8) 

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p034r033/images'; % (8) 

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p036r038/images'; % (8) 

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p022r033/images'; % (8) 

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p039r026/images'; % (8) 

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p013r029/images'; % (8) 

% v_input.l_dir='/projectnb/landsat/projects/LCMAP/p031r027/images'; % (8) 
v_input.l_dir = pwd % (6) 

% 2. ground truth time interval for continuous classification
% v_input.gt=[datenum(2005,1,0),datenum(2007,12,31)]; % Boston
% v_input.gt=[datenum(2001,4,17),datenum(2001,10,26)]; % Amazon
% v_input.gt=[datenum(2000,1,1),datenum(2000,12,31)]; % LCMS
% v_input.gt=[datenum(2000,1,1),datenum(2000,12,31)]; % Composite yearly
% v_input.gt=[datenum(2000,6,1),datenum(2000,9,1)]; % Composite seasonly (June~August)
% v_input.gt=[datenum(2008,1,1),datenum(2009,12,31)]; % Guangzhou
% v_input.gt = [datenum(1999,1,1),datenum(2001,12,31)]; % Trends
% v_input.gt = [datenum(2006,1,1),datenum(2006,12,31)]; % NLCD 2006


%% Constants:
% number of maximum coefficients
v_input.num_c = 8; 
% name of ground truth land cover map
v_input.name_roi = 'example_img'; 
% v_input.name_roi = 'training_data'; 
% folder name of all CCDC resultls 
v_input.name_rst = 'TSFitMap';
% folder name of CCDC maps
v_input.name_map = 'CCDCMap';
% folder name of CCDC predicted surf ref
v_input.name_pre = 'PredictAll';
% v_input.name_roi='class_map'; 
% get to the data folder
% cd(v_input.l_dir);
% read in image dimension from reference image
% [jiDim,jiul,resolu,zc] = envihdrread(v_input.name_roi);

% get num of total folders start with "L"
imf=dir('L*'); % folder names
% filter for Landsat folders
imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
imf = [imf{:}];
imf = vertcat(imf{:});
% sort according to yeardoy
yeardoy = str2num(imf(:, 10:16));
[~, sort_order] = sort(yeardoy);
imf = imf(sort_order, :);
% name of the first stacked image
filename = dir([imf(1,:),'/','L*MTLstack']); 
% read in dimension and zone number of the data
[jiDim,jiul,resolu,zc] = envihdrread([imf(1,:),'/',filename.name]);

% jiDim = [1000 1000];
% jiul = [0 0];
% resolu = 30;
% zc = 0;
% dimension of image [row,col]
v_input.ijdim = [jiDim(2),jiDim(1)];
% upper left coordinates
v_input.jiul = jiul;
% resolution of the image
v_input.resolu = resolu;
% zone code
v_input.zc = zc;
% number of bands (7 spectral bands + Fmask)
v_input.nbands = 8;
end