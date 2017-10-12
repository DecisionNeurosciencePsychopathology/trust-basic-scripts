function s = censorComputer(id,s,block)
%Censor entire blocks if max movement is greater than 5 per block
%Load in the data
path_of_conditions_table = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/OrderTrustee.mat';
load(path_of_conditions_table);
%load('fd_max.mat') %Change path as needed

if ischar(id)
    id = str2double(id);
end

%First get index of id or kick out
id_idx=find(ismember(OrderTrustee.ID,id), 1);
if isempty(id_idx)
    warning('Subject is not on the block list')
    return
end

censor_limit = 4;
row_name = ['Block' num2str(block)];

%Censor the entire block if greater than the censor limit
if OrderTrustee.(row_name)(id_idx)>=censor_limit
    s.(['regressors' num2str(block)]).to_censor = ...
        logical(s.(['regressors' num2str(block)]).to_censor .* 0);
    
end

