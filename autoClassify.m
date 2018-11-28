function autoClassify(task,ntasks,varargin)
%AUTOCLASSIFICATION runs for Continous Classification.
% 
%AUTOCLASSIFY(1,1) automatically goes continuous classification at 
%                        the current directory using 1 core.
%
% Specific parameters
% ------------------------
%   'task'              CPU identification number for parallel computing. *REQUIRED.
%   'ntasks'            Total number of CPU used. *REQUIRED.
%   'NTrees'            Number of trees grown. Default is 500.
%   'CCDCDir'           Directory of CCDC change detection results.  
%                       Default is the path to current folder.
%
% History
% ------------------------
% Revisions: $ Date: 11/25/2015 $ Copyright: Zhe Zhu
% Version 1.1: Combine classification and change detection output (07/05/2015)
% Version 1.0: Fast classification for each line (01/10/2015)
%
% Author: Zhe Zhu (zhe.zhu#ttu.edu)
%         Shi Qiu (shi.qiu#ttu.edu)
% Date: 28. Jun, 2018

%     task = 1;
%     ntasks = 1;
    
      %% get parameters from inputs
    % defaults
    % where the all CCDC change detection results are
    dir_cur = pwd;
 
    % number of trees grown.
    ntrees = 500;
    % request the user's inputs
    p = inputParser;
    p.FunctionName = 'classifyParas';
    % optional
    % default values.
    addParameter(p,'CCDCDir',dir_cur);
    addParameter(p,'NTrees',ntrees);
    
    parse(p,varargin{:});
    ccdc_dir = p.Results.CCDCDir;
    ntrees = p.Results.NTrees; % number of trees grown.
    
    
    fprintf('Start...\n');
    %
    % get image parameters automatically
    imf=dir(fullfile(ccdc_dir,'L*')); % folder names
%             imf=dir('L*'); % folder names
    [nrows,~,nbands] = autoPara(imf);

    %% Constants: 
    % number of trees
% %     ntrees = 500;
    % version of CCDC
    ccdc_v = 1.1;
    % number of coefficient
    num_c = 8;

    % log for CCDC Change paramters and versions
    % report only for the first task
    if task == 1
        fileID = fopen(fullfile(ccdc_dir,'CCDC_Classification_log.txt'),'w');
        % write location of image stack
        fprintf(fileID,'Image location = %s\r\n',ccdc_dir);
        % write number of images used
        fprintf(fileID,'Number of trees = %d\r\n',ntrees);
        % CCDC Version
        fprintf(fileID,'CCDC Classification Version = %.2f\r\n',ccdc_v);
        % updates
        fprintf(fileID,'******************************************************************************************************\r\n');
        fprintf(fileID,'Revisions: $ Date: 11/20/2015 $ Copyright: Zhe Zhu\r\n');
        fprintf(fileID,'Version 1.1   Combine classification and change detection output (07/05/2015)\r\n');
        fprintf(fileID,'Version 1.0   Fast classification for each line (01/10/2015)\r\n');
        fprintf(fileID,'******************************************************************************************************\r\n');
        fclose(fileID);
    end

    % load in RF models 
    load( fullfile(ccdc_dir,'modelRF')); % change name if neccessary

    % make TSFitMap folder for storing coefficients
    n_result = 'TSFitMap';
    if isempty(dir(fullfile(ccdc_dir,n_result)))
        mkdir(dir(fullfile(ccdc_dir,n_result)));
    end

    % prepare the irows for task for ALL rows
    irows = zeros(1,1);
    i = 0;
    while task + ntasks*i <= nrows
       irows(i+1) = task + ntasks*i;
       i = i+1;
    end

    for i = 1:length(irows)   
        try
            % if successfully loaded => skip
%                     load([dir_l,'/',n_result,'/','record_change',sprintf('%d',irows(i)),'.mat']);
            load(fullfile(ccdc_dir,n_result,['record_change',sprintf('%d',irows(i)),'.mat']));
        catch me
            % not exist or corrupt
%                     fprintf('Missing the %dth row!\n',irows(i));

            fprintf('Missing the %dth row!\n',irows(i));

            continue
        end

        fprintf('Processing the %dth row.\n',irows(i));

%                 fprintf('Processing the %dth row\n',irows(i));
        % Continous Classfication Done for a Line of timeseries pixels
        Class_Line1_1(ccdc_dir,n_result,irows(i),modelRF,num_c,nbands,ntrees);
    end
    fprintf('Done!\n');
end

