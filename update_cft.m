function update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c)

% determine the time series model
if i_span < mid_num_c*n_times
    % start with 4 coefficients model
    update_num_c = min(min_num_c,num_c);
elseif i_span < max_num_c*n_times
    % start with 6 coefficients model
    update_num_c = min(mid_num_c,num_c);
else
    % start with 8 coefficients model
    update_num_c =  min(max_num_c,num_c);
end

end