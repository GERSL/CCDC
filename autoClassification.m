function autoClassification(task,ntasks,varargin)
%AUTOCLASSIFICATION runs for Continous Classification.
% 
%AUTOCLASSIFICATION(1,1) automatically goes continuous classification at 
%                        the current directory using 1 core.
%
% Specific parameters
% ------------------------
%   'task'              CPU identification number for parallel computing. *REQUIRED.
%   'ntasks'            Total number of CPU used. *REQUIRED.
%   'NTrees'            Number of trees grown. Default is 500.
%   'CCDCDirectory'     Directory of CCDC change detection results.  
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
    addParameter(p,'CCDCDirectory',dir_cur);
    addParameter(p,'NTrees',ntrees);
    
    parse(p,varargin{:});
    ccdc_dir = p.Results.CCDCDirectory;
    ntrees = p.Results.NTrees; % number of trees grown.
    
    
    fprintf('Start...');
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

            fprintf('Missing the %dth row!',irows(i));

            continue
        end

        fprintf('Processing the %dth row.',irows(i));

%                 fprintf('Processing the %dth row\n',irows(i));
        % Continous Classfication Done for a Line of timeseries pixels
        Class_Line1_1(ccdc_dir,n_result,irows(i),modelRF,num_c,nbands,ntrees);
    end
    fprintf('Done!\n');
end

function Class_Line1_1(dir_l,n_rst,nrow,model,num_c,nbands,ntrees)
% This is the classification algorithm for classifying all map pixels by lines
%
% Revisions: $ Date: 11/25/2015 $ Copyright: Zhe Zhu
% Version 1.1: Combine classification and change detection output (07/05/2015)
% Version 1.0: Fast classification for each line (01/10/2015)

% fprintf('Processing the %d row\n',i);
% load([dir_l,'/',n_rst,'/','record_change',num2str(nrow)]);
load(fullfile(dir_l,n_rst,['record_change',num2str(nrow)]));

% do not classify if class exist
% if isfield(rec_cg,'class') == 0
    
    % total number of time series models (including empty ones)
    num_ts =  length(rec_cg); %#ok<*NODEF>
    % initiate Ids
    IDs_rec = zeros(1,num_ts);
    
    for i_model = 1:num_ts
        % record the time series ids that have recorded values
        % for adding classification results
        if ~isempty(rec_cg(i_model).pos)
            IDs_rec(i_model) = 1;
        end
    end
    
    % start time
    t_start = [rec_cg.t_start];
    % end time
    t_end = [rec_cg.t_end];
    % position
    pos = [rec_cg.pos];
    % rmse
    rmse = [rec_cg.rmse]; % each curve has nbands rmse
    % number of valid curves per line
    num_s = sum(IDs_rec);
    
    if num_s > 0 % has more than one curves exist for each line
        % model coefficients
        tmp = [rec_cg.coefs];
        % prepare for classification inputs
        Xclass = zeros(num_s,(num_c+1)*(nbands-1));
        
        for icol = 1:num_s;
            % coefficients from the 7 bands
            i_tmp = tmp(:,((icol-1)*(nbands-1)+1):(nbands-1)*icol);
            % modified constant as inputs
            i_tmp(1,:) = i_tmp(1,:)+(t_start(icol)+t_end(icol))*i_tmp(2,:)/2;
            % input ready!
            Xclass(icol,:) = [reshape(rmse(((icol-1)*(nbands-1)+1):(nbands-1)*icol),nbands-1,1);i_tmp(:)];
        end
        % classify the whole line
        [map,votes] = classRF_predict(Xclass,model,ntrees); % class
    end
    
    % add a new component "class" to rec_cg
    if sum(IDs_rec == 0) > 0
        IDs_add = find(IDs_rec == 0);
        for i = 1:length(IDs_add)
            rec_cg(IDs_add(i)).class = [];
            rec_cg(IDs_add(i)).classQA = [];
        end
    end
    
    if sum(IDs_rec == 1) > 0
        IDs_add = find(IDs_rec == 1);
        for i = 1:length(IDs_add)
            rec_cg(IDs_add(i)).class = map(i);
            
            % largest number of votes
            [max_v1,max_id] = max(votes(i,:));
            % make this smallest
            votes(i,max_id) = 0;
            % second largest number of votes
            max_v2 = max(votes(i,:));
            % provide unsupervised ensemble margin as QA
            rec_cg(IDs_add(i)).classQA = 100*(max_v1-max_v2)/ntrees;
        end
    end
    
    % updated "record_change*.mat"
%     save([dir_l,'/',n_rst,'/','record_change',num2str(nrow)],'rec_cg');
    save(fullfile(dir_l,n_rst,['record_change',num2str(nrow)]),'rec_cg');
% end

end % end of the function
