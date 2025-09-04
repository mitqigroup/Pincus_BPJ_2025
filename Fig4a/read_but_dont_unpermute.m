function data_structure = read_but_dont_unpermute(data_path, num_tasks)

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
        'parameter_set', 'phi_end_nonlinear');
else
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
        'parameter_set', 'phi_end_nonlinear_self');
%         load(load_name);
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
    
    parameter_set = data.parameter_set;
    save(fullfile(data_path,'combined.mat'), ...
        'E_all_toroid', 'alpha_A_vals_toroid', 'alpha_B_vals_toroid',...
        'S_A_vals_toroid', 'S_B_vals_toroid', 'Sigma_vals_toroid',...
        'rho_vals', 'E_all_linear','alpha_A_vals_linear','alpha_B_vals_linear',...
        'S_A_vals_linear','S_B_vals_linear','Sigma_vals_linear',...
        'E_all_nonlinear','alpha_A_vals_nonlinear','alpha_B_vals_nonlinear',...
        'S_A_vals_nonlinear','S_B_vals_nonlinear','Sigma_vals_nonlinear',...
        'parameter_set', 'phi_end_nonlinear');
end
toc

data_structure.parameter_set = parameter_set;
data_structure.E_all_toroid = E_all_toroid;
data_structure.alpha_A_vals_toroid = alpha_A_vals_toroid;
data_structure.alpha_B_vals_toroid = alpha_B_vals_toroid;
data_structure.S_A_vals_toroid = S_A_vals_toroid;
data_structure.S_B_vals_toroid = S_B_vals_toroid;
data_structure.Sigma_vals_toroid = Sigma_vals_toroid;
data_structure.rho_vals = rho_vals;

data_structure.E_all_linear = E_all_linear;
data_structure.alpha_A_vals_linear = alpha_A_vals_linear;
data_structure.alpha_B_vals_linear = alpha_B_vals_linear;
data_structure.S_A_vals_linear = S_A_vals_linear;
data_structure.S_B_vals_linear = S_B_vals_linear;
data_structure.Sigma_vals_linear = Sigma_vals_linear;

data_structure.E_all_nonlinear = E_all_nonlinear;
data_structure.alpha_A_vals_nonlinear = alpha_A_vals_nonlinear;
data_structure.alpha_B_vals_nonlinear = alpha_B_vals_nonlinear;
data_structure.S_A_vals_nonlinear = S_A_vals_nonlinear;
data_structure.S_B_vals_nonlinear = S_B_vals_nonlinear;
data_structure.Sigma_vals_nonlinear = Sigma_vals_nonlinear;
if isfield(data, 'phi_end_nonlinear_self')
    data_structure.phi_end_nonlinear = phi_end_nonlinear;
end

end