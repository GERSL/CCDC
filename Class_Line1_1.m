
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