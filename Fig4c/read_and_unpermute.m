function data_structure = read_and_unpermute(data_path, num_tasks, varargin)

% I should probably re-write this programattically using dynamic field
% names, basically vars_list = ['E_all_toroid', 'Sigma_vals_toroid', ...]
% and then data_structure.(vars_list(ii)) = vals, with different types to
% unpermute them differently. That way it will be easier to extend this to
% adding extra physical parameters.

% parse inputs and set defaults
args = varargin;
nargs = numel(args);
sort_vals = false;
only_nonlinear = false;
k = 1;
while k<=nargs
    if strcmpi(args{k},'SortValues')||strcmpi(args{k},'Sort_values')
        sort_vals = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'OnlyNonlinear')||strcmpi(args{k},'Only_Nonlinear')
        only_nonlinear = args{k+1};
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
        'E_all_toroid', 'alpha_A_vals_toroid', 'alpha_B_vals_toroid',...
        'S_A_vals_toroid', 'S_B_vals_toroid', 'Sigma_vals_toroid',...
        'rho_vals', 'E_all_linear','alpha_A_vals_linear','alpha_B_vals_linear',...
        'S_A_vals_linear','S_B_vals_linear','Sigma_vals_linear',...
        'E_all_nonlinear','alpha_A_vals_nonlinear','alpha_B_vals_nonlinear',...
        'S_A_vals_nonlinear','S_B_vals_nonlinear','Sigma_vals_nonlinear',...
        'parameter_set', 'phi_end_nonlinear', 'phi_vals');
else
    if ~only_nonlinear
        E_all_toroid = [];
        alpha_A_vals_toroid = [];
        alpha_B_vals_toroid = [];
        S_A_vals_toroid = [];
        S_B_vals_toroid = [];
        Sigma_vals_toroid = [];
        rho_vals = [];
        
        E_all_linear = [];
        alpha_A_vals_linear = [];
        alpha_B_vals_linear = [];
        S_A_vals_linear = [];
        S_B_vals_linear = [];
        Sigma_vals_linear = [];
    end
    
    E_all_nonlinear = [];
    alpha_A_vals_nonlinear = [];
    alpha_B_vals_nonlinear = [];
    S_A_vals_nonlinear = [];
    S_B_vals_nonlinear = [];
    Sigma_vals_nonlinear = [];
    phi_end_nonlinear = [];

    for ii=0:num_tasks-1
        load_name = sprintf('data/task%i_results.mat', ii )
        data = load(fullfile(data_path,load_name), ...
        'E_all_toroid_self', 'alpha_A_vals_toroid_self', 'alpha_B_vals_toroid_self',...
        'S_A_vals_toroid_self', 'S_B_vals_toroid_self', 'Sigma_vals_toroid_self',...
        'rho_vals_self', 'E_all_linear_self','alpha_A_vals_linear_self','alpha_B_vals_linear_self',...
        'S_A_vals_linear_self','S_B_vals_linear_self','Sigma_vals_linear_self',...
        'E_all_nonlinear_self','alpha_A_vals_nonlinear_self','alpha_B_vals_nonlinear_self',...
        'S_A_vals_nonlinear_self','S_B_vals_nonlinear_self','Sigma_vals_nonlinear_self',...
        'parameter_set', 'phi_end_nonlinear_self', 'phi_vals');
%         load(load_name);
        if ~only_nonlinear
            E_all_toroid = cat(3,E_all_toroid, data.E_all_toroid_self);
            alpha_A_vals_toroid = horzcat(alpha_A_vals_toroid, data.alpha_A_vals_toroid_self);
            alpha_B_vals_toroid = horzcat(alpha_B_vals_toroid, data.alpha_B_vals_toroid_self);
            S_A_vals_toroid = horzcat(S_A_vals_toroid, data.S_A_vals_toroid_self);
            S_B_vals_toroid = horzcat(S_B_vals_toroid, data.S_B_vals_toroid_self);
            Sigma_vals_toroid = horzcat(Sigma_vals_toroid, data.Sigma_vals_toroid_self);
            rho_vals = horzcat(rho_vals, data.rho_vals_self);
            
            E_all_linear = cat(3,E_all_linear, data.E_all_linear_self);
            alpha_A_vals_linear = horzcat(alpha_A_vals_linear, data.alpha_A_vals_linear_self);
            alpha_B_vals_linear = horzcat(alpha_B_vals_linear, data.alpha_B_vals_linear_self);
            S_A_vals_linear = horzcat(S_A_vals_linear, data.S_A_vals_linear_self);
            S_B_vals_linear = horzcat(S_B_vals_linear, data.S_B_vals_linear_self);
            Sigma_vals_linear = horzcat(Sigma_vals_linear , data.Sigma_vals_linear_self);
        end
        
        E_all_nonlinear = cat(3,E_all_nonlinear, data.E_all_nonlinear_self);
        alpha_A_vals_nonlinear = horzcat(alpha_A_vals_nonlinear, data.alpha_A_vals_nonlinear_self);
        alpha_B_vals_nonlinear = horzcat(alpha_B_vals_nonlinear, data.alpha_B_vals_nonlinear_self);
        S_A_vals_nonlinear = horzcat(S_A_vals_nonlinear, data.S_A_vals_nonlinear_self);
        S_B_vals_nonlinear = horzcat(S_B_vals_nonlinear, data.S_B_vals_nonlinear_self);
        Sigma_vals_nonlinear = horzcat(Sigma_vals_nonlinear, data.Sigma_vals_nonlinear_self);
        if isfield(data, 'phi_end_nonlinear_self')
            phi_end_nonlinear = horzcat(phi_end_nonlinear, data.phi_end_nonlinear_self);
        end

    end

    phi_vals = data.phi_vals;
    
    parameter_set = data.parameter_set;
    if only_nonlinear
        save(fullfile(data_path,'combined.mat'), ...
            'E_all_nonlinear','alpha_A_vals_nonlinear','alpha_B_vals_nonlinear',...
            'S_A_vals_nonlinear','S_B_vals_nonlinear','Sigma_vals_nonlinear',...
            'parameter_set', 'phi_end_nonlinear', 'phi_vals');
    else
        save(fullfile(data_path,'combined.mat'), ...
            'E_all_toroid', 'alpha_A_vals_toroid', 'alpha_B_vals_toroid',...
            'S_A_vals_toroid', 'S_B_vals_toroid', 'Sigma_vals_toroid',...
            'rho_vals', 'E_all_linear','alpha_A_vals_linear','alpha_B_vals_linear',...
            'S_A_vals_linear','S_B_vals_linear','Sigma_vals_linear',...
            'E_all_nonlinear','alpha_A_vals_nonlinear','alpha_B_vals_nonlinear',...
            'S_A_vals_nonlinear','S_B_vals_nonlinear','Sigma_vals_nonlinear',...
            'parameter_set', 'phi_end_nonlinear', 'phi_vals');
    end

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
    if ~only_nonlinear
        parameter_set_unpermuted(kk,:) = parameter_set(jj,:);
        E_all_toroid_unpermuted(:,:,kk) = E_all_toroid(:,:,jj);
        alpha_A_vals_toroid_unpermuted(:,kk) = alpha_A_vals_toroid(:,jj);
        alpha_B_vals_toroid_unpermuted(:,kk) = alpha_B_vals_toroid(:,jj);
        S_A_vals_toroid_unpermuted(:,kk) = S_A_vals_toroid(:,jj);
        S_B_vals_toroid_unpermuted(:,kk) = S_B_vals_toroid(:,jj);
        Sigma_vals_toroid_unpermuted(:,kk) = Sigma_vals_toroid(:,jj);
        rho_vals_unpermuted(:,kk) = rho_vals(:,jj);
        
        E_all_linear_unpermuted(:,:,kk) = E_all_linear(:,:,jj);
        alpha_A_vals_linear_unpermuted(:,kk) = alpha_A_vals_linear(:,jj);
        alpha_B_vals_linear_unpermuted(:,kk) = alpha_B_vals_linear(:,jj);
        S_A_vals_linear_unpermuted(:,kk) = S_A_vals_linear(:,jj);
        S_B_vals_linear_unpermuted(:,kk) = S_B_vals_linear(:,jj);
        Sigma_vals_linear_unpermuted(:,kk) = Sigma_vals_linear(:,jj);
    end
    
    E_all_nonlinear_unpermuted(:,:,kk) = E_all_nonlinear(:,:,jj);
    alpha_A_vals_nonlinear_unpermuted(:,kk) = alpha_A_vals_nonlinear(:,jj);
    alpha_B_vals_nonlinear_unpermuted(:,kk) = alpha_B_vals_nonlinear(:,jj);
    S_A_vals_nonlinear_unpermuted(:,kk) = S_A_vals_nonlinear(:,jj);
    S_B_vals_nonlinear_unpermuted(:,kk) = S_B_vals_nonlinear(:,jj);
    Sigma_vals_nonlinear_unpermuted(:,kk) = Sigma_vals_nonlinear(:,jj);
%     if isfield(data, 'phi_end_nonlinear_self')
    if exist('phi_end_nonlinear', 'var')&&~isempty(phi_end_nonlinear)
        phi_end_nonlinear_unpermuted(:,kk) = phi_end_nonlinear(:,jj);
    end
end

if sort_vals
    [phi_vals,sortID] = sort(phi_vals,'ascend');
    if ~only_nonlinear
        E_all_toroid_unpermuted = E_all_toroid_unpermuted(:,sortID,:);
        alpha_A_vals_toroid_unpermuted = alpha_A_vals_toroid_unpermuted(sortID,:);
        alpha_B_vals_toroid_unpermuted = alpha_B_vals_toroid_unpermuted(sortID,:);
        S_A_vals_toroid_unpermuted = S_A_vals_toroid_unpermuted(sortID,:);
        S_B_vals_toroid_unpermuted = S_B_vals_toroid_unpermuted(sortID,:);
        Sigma_vals_toroid_unpermuted = Sigma_vals_toroid_unpermuted(sortID,:);
        rho_vals_unpermuted = rho_vals_unpermuted(sortID,:);
        
        E_all_linear_unpermuted = E_all_linear_unpermuted(:,sortID,:);
        alpha_A_vals_linear_unpermuted = alpha_A_vals_linear_unpermuted(sortID,:);
        alpha_B_vals_linear_unpermuted = alpha_B_vals_linear_unpermuted(sortID,:);
        S_A_vals_linear_unpermuted = S_A_vals_linear_unpermuted(sortID,:);
        S_B_vals_linear_unpermuted = S_B_vals_linear_unpermuted(sortID,:);
        Sigma_vals_linear_unpermuted = Sigma_vals_linear_unpermuted(sortID,:);
    end
    
    E_all_nonlinear_unpermuted = E_all_nonlinear_unpermuted(:,sortID,:);
    alpha_A_vals_nonlinear_unpermuted = alpha_A_vals_nonlinear_unpermuted(sortID,:);
    alpha_B_vals_nonlinear_unpermuted = alpha_B_vals_nonlinear_unpermuted(sortID,:);
    S_A_vals_nonlinear_unpermuted = S_A_vals_nonlinear_unpermuted(sortID,:);
    S_B_vals_nonlinear_unpermuted = S_B_vals_nonlinear_unpermuted(sortID,:);
    Sigma_vals_nonlinear_unpermuted = Sigma_vals_nonlinear_unpermuted(sortID,:);
%     if isfield(data, 'phi_end_nonlinear_self')
    if exist('phi_end_nonlinear', 'var')&&~isempty(phi_end_nonlinear)
        phi_end_nonlinear_unpermuted = phi_end_nonlinear_unpermuted(sortID,:);
    end
end

if ~only_nonlinear
    data_structure.parameter_set = parameter_set_unpermuted;
    data_structure.E_all_toroid = E_all_toroid_unpermuted;
    data_structure.alpha_A_vals_toroid = alpha_A_vals_toroid_unpermuted;
    data_structure.alpha_B_vals_toroid = alpha_B_vals_toroid_unpermuted;
    data_structure.S_A_vals_toroid = S_A_vals_toroid_unpermuted;
    data_structure.S_B_vals_toroid = S_B_vals_toroid_unpermuted;
    data_structure.Sigma_vals_toroid = Sigma_vals_toroid_unpermuted;
    data_structure.rho_vals = rho_vals_unpermuted;
    
    data_structure.E_all_linear = E_all_linear_unpermuted;
    data_structure.alpha_A_vals_linear = alpha_A_vals_linear_unpermuted;
    data_structure.alpha_B_vals_linear = alpha_B_vals_linear_unpermuted;
    data_structure.S_A_vals_linear = S_A_vals_linear_unpermuted;
    data_structure.S_B_vals_linear = S_B_vals_linear_unpermuted;
    data_structure.Sigma_vals_linear = Sigma_vals_linear_unpermuted;
end

data_structure.E_all_nonlinear = E_all_nonlinear_unpermuted;
data_structure.alpha_A_vals_nonlinear = alpha_A_vals_nonlinear_unpermuted;
data_structure.alpha_B_vals_nonlinear = alpha_B_vals_nonlinear_unpermuted;
data_structure.S_A_vals_nonlinear = S_A_vals_nonlinear_unpermuted;
data_structure.S_B_vals_nonlinear = S_B_vals_nonlinear_unpermuted;
data_structure.Sigma_vals_nonlinear = Sigma_vals_nonlinear_unpermuted;
if exist('phi_end_nonlinear', 'var')&&~isempty(phi_end_nonlinear)
    data_structure.phi_end_nonlinear = phi_end_nonlinear_unpermuted;
end
data_structure.phi_vals = phi_vals;
data_structure.is_sorted = sort_vals;

end