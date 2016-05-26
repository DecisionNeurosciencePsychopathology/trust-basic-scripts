function s = censorMovement(id,s,block)
%Censor entire blocks if max movement is greater than 5 per block
%Load in the data
load('fd_max.mat') %Change path as needed

%First get index of id or kick out
id_idx=find(ismember(fd_max.Subjects,id), 1);
if isempty(id_idx)
    warning('Subject is not on censor block list')
    return
end

censor_limit = 5;
row_name = ['Max' num2str(block)];

%Censor the entire block if greater than the censor limit
if fd_max.(row_name)(id_idx)>=censor_limit
    s.(['regressors' num2str(block)]).to_censor = ...
        logical(s.(['regressors' num2str(block)]).to_censor .* 0);
    
end

