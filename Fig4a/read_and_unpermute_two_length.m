function data_structure = read_and_unpermute_two_length(data_path, num_tasks, varargin)

% parse inputs and set defaults
args = varargin;
nargs = numel(args);
sort_vals = false;
k = 1;
while k<=nargs
    if strcmpi(args{k},'SortValues')||strcmpi(args{k},'Sort_values')
        sort_vals = args{k+1};
        k = k+1;
    else
        error('intput not recognised!');
    end
    k = k+1;
end

data = struct;

tic
if isfile(fullfile(data_path,'combined.mat'))
    load(fullfile(data_path,'combined.mat'), ...
        'E', 'aA', 'aB', 'phi_vals', 'parameter_set');
else
    E = [];
    aA = [];
    aB = [];

    for ii=0:num_tasks-1
        load_name = sprintf('data/task%i_results.mat', ii )
        data = load(fullfile(data_path,load_name), ...
        'E_self', 'aA_self', 'aB_self', 'phi_vals', 'parameter_set');
%         load(load_name);
        E = vertcat(E, data.E_self);
        aA = vertcat(aA, data.aA_self);
        aB = vertcat(aB, data.aB_self);

    end

    phi_vals = data.phi_vals;
    
    parameter_set = data.parameter_set;
    save(fullfile(data_path,'combined.mat'), ...
        'E', 'aA', 'aB', 'phi_vals', 'parameter_set');
end
toc

% randomly shuffle the parameter set so that we have a more even distribution between
% each core on our node. We can unshuffle it later on if we want to. We MUST set the
% random number generator to the same value, since it must be identical between cores
rng('default');
rng(1);

size_params = size(parameter_set);

permutation_array = randperm(size_params(1));

% parameter_set_permuted = parameter_set;
% % parameter_set_unpermuted = parameter_set;
% 
% for ii=1:size_params(2)
%     parameter_set_permuted(:,ii) = parameter_set(permutation_array,ii);
% end
% parameter_set_original = parameter_set;
% parameter_set = parameter_set_permuted;

for jj = 1:length(permutation_array)
    kk = permutation_array(jj);
    parameter_set_unpermuted(kk,:) = parameter_set(jj,:);
    E_unpermuted(kk,:) = E(jj,:);
    aA_unpermuted(kk,:) = aA(jj,:);
    aB_unpermuted(kk,:) = aB(jj,:);
end

if sort_vals
    [phi_vals,sortID] = sort(phi_vals,'ascend');
    E_unpermuted = E_unpermuted(:,sortID);
    aA_unpermuted = aA_unpermuted(:,sortID);
    aB_unpermuted = aB_unpermuted(:,sortID);
end

data_structure.parameter_set = parameter_set_unpermuted;
data_structure.E = E_unpermuted;
data_structure.aA = aA_unpermuted;
data_structure.aB = aB_unpermuted;
data_structure.phi_vals = phi_vals;
data_structure.is_sorted = sort_vals;

end