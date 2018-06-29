function [nrows,ncols,nbands,jiUL,resolu,zc,num_imgs] = autoPara(imf)
% Matlab code for read paramters of files automatically
% Inputs:
% imf:   folder names from dir

% number of images
num_imgs = size(imf,1);
% filter for Landsat folders
imf_path = imf.folder;
imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
imf = [imf{:}];
imf = vertcat(imf{:});
% name of the first stacked image
filename = dir(fullfile(imf_path(1,:),imf(1,:),'L*stack')); 
% read in dimension and zone number of the data
[jiDim,jiUL,resolu,zc,nbands] = envihdrread(fullfile(imf_path(1,:),imf(1,:),filename.name));
% number of nrows processed
nrows = jiDim(2);
% number of pixels procesed per line
ncols = jiDim(1); 
end