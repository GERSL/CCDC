%**************************************************************
%* mex interface to Andy Liaw et al.'s C code (used in R package randomForest)
%* Added by Abhishek Jaiantilal ( abhishek.jaiantilal@colorado.edu )
%* License: GPLv2
%* Version: 0.02
%
% Calls Classification Random Forest
% A wrapper matlab file that calls the mex file
% This does prediction given the data and the model file
% Options depicted in predict function in http://cran.r-project.org/web/packages/randomForest/randomForest.pdf
%**************************************************************
%function [Y_hat votes] = classRF_predict(X,model, extra_options)
% requires 2 arguments
% X: data matrix
% model: generated via classRF_train function
% extra_options.predict_all = predict_all if set will send all the prediction. 
%
%
% Returns
% Y_hat - prediction for the data
% votes - unnormalized weights for the model
% prediction_per_tree - per tree prediction. the returned object .
%           If predict.all=TRUE, then the individual component of the returned object is a character
%           matrix where each column contains the predicted class by a tree in the forest.
% proximity_ts - proximity of training to the test set
% nodes - Should the terminal node indicators (an n by ntree matrix) be return? If so, it is
%         in the “nodes” attribute of the returned object.
% 
% Not yet implemented
% proximity

function [Y_new, votes, prediction_per_tree,proximity_ts,nodes] = classRF_predict(X,model, extra_options)
    
    if nargin<2
		error('need atleast 2 parameters,X matrix and model');
    end
    
    if exist('extra_options','var')
        if isfield(extra_options,'predict_all') 
            predict_all = extra_options.predict_all;
        end
        if isfield(extra_options,'proximity') 
            proximity = extra_options.proximity;
        end
        if isfield(extra_options,'nodes') 
            nodes = extra_options.nodes;
        end
    end
    
    if ~exist('predict_all','var'); predict_all=0;end
            
    if ~exist('proximity','var'); proximity=0;end
        
    if ~exist('nodes','var'); nodes=0;end
    
    if isfield(model, 'categorical_feature')
        % have to map prediction array to the correct categories
        for i=1:size(X,2)
            if model.categorical_feature(i) 
                tmp_uniques_in_feature = model.orig_uniques_in_feature{i};
                tmp_mapped_uniques_in_feature = model.mapped_uniques_in_feature{i};
                X_loc = X(:,i); %cannot change the original array which may cause chained change of categories to something totally wrong
                for j=1:length(tmp_uniques_in_feature)
                    indices_to_change = find( X(:,i) == tmp_uniques_in_feature(j) );
                    X_loc(indices_to_change) = tmp_mapped_uniques_in_feature(j);
                end
                X(:,i) = X_loc;
            end
        end
        ncat = model.ncat;
    else
        ncat = ones(1,size(X,2));
    end    
    
    maxcat = max(ncat);
	[Y_hat,prediction_per_tree,votes,proximity_ts,nodes] = mexClassRF_predict(X',model.nrnodes,model.ntree,model.xbestsplit,model.classwt,model.cutoff,model.treemap,model.nodestatus,model.nodeclass,model.bestvar,model.ndbigtree,model.nclass, predict_all, proximity, nodes, int32(ncat), maxcat );
	%keyboard
    votes = votes';
    
    clear mexClassRF_predict
    
    Y_new = double(Y_hat);
    new_labels = model.new_labels;
    orig_labels = model.orig_labels;
    
    for i=1:length(orig_labels)
        Y_new(find(Y_hat==new_labels(i)))=Inf;
        Y_new(isinf(Y_new))=orig_labels(i);
    end
    
    1;
    