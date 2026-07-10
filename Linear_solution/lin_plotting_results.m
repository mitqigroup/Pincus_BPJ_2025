clear variables

tic
if isfile('combined.mat')
    load('combined.mat')
else
    num_tasks = 48;
    solution_set_total = {};
    E_all_total = [];
    alpha_A_vals_total = [];
    alpha_B_vals_total = [];
    phi_vals_total = [];
    h_phi_vals_total = [];
    S_A_vals_total = [];
    S_B_vals_total = [];
    Sigma_vals_total = [];
    for ii=0:num_tasks-1
        load_name = sprintf('data/task%i_results.mat', ii )
        data = load(load_name);
        solution_set_total = [solution_set_total, data.solution_set];
        E_all_total = cat(2,E_all_total, data.E_all);
        alpha_A_vals_total = horzcat(alpha_A_vals_total, data.alpha_A_vals);
        alpha_B_vals_total = horzcat(alpha_B_vals_total, data.alpha_B_vals);
        phi_vals_total = horzcat(phi_vals_total, data.phi_vals);
        h_phi_vals_total = horzcat(h_phi_vals_total, data.h_phi_vals);
        S_A_vals_total = horzcat(S_A_vals_total, data.S_A_vals);
        S_B_vals_total = horzcat(S_B_vals_total, data.S_B_vals);
        Sigma_vals_total = horzcat(Sigma_vals_total, data.Sigma_vals);
    end
    
    parameter_set = data.parameter_set;
    save('combined.mat')
end
toc

%% base case
% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = [0.05, 0.3, 1];                 	   % um
sigma_vals = logspace(-3,-1,24);               % surface fraction
%d_vals = sqrt(R_vals.^2./sigma_vals);         % um
kD_vals = 300/10^12*1e9;                       % picoJ/um^2
kD_base = 300/10^12*1e9;                       % picoJ/um^2
zeta_vals = [0.001,0.005,0.02];                % dimensionless
epsilon_vals = -zeta_vals*kD_base;             % picoJ/um^2
n0_vals = 1;                                   % fraction
kappa_vals = [1e-19*1e12, 1e-18*1e12];         % picoJ
% kappa_vals = logspace(-21,-15, 60)*1e12;
alpha_i_vals = [0,0.013];                  % fraction
% other constants
% N = 3e5;

ii=0;
slice = [];
for rr = 1:length(R_vals)
    for ss = 1:length(sigma_vals)
        for kk = 1:length(kD_vals)
            for ee = 1:length(epsilon_vals)
                for nn = 1:length(n0_vals)
                    for pp = 1:length(kappa_vals)
                        for aa = 1:length(alpha_i_vals)
                            ii = ii+1;
                            if ee==3&&rr==2&&pp==1&&aa==1
                                slice = [slice, ii];
                            end
                        end
                    end
                end
            end
        end
    end
end
size_params = size(parameter_set);
base_index = slice(1);
% base_index = 1;
% step = 1;
% % if using R
% slice = base_index:step:size_params(1);
% if using alpha_i
% slice = (base_index-1)*step+1:base_index*step;
% base_index = (base_index-1)*step+1;
epsilon = parameter_set(base_index,1);
n0 = parameter_set(base_index,2);
d = parameter_set(base_index,3);
R = parameter_set(base_index,4);
kD = parameter_set(base_index,5);
kappa = parameter_set(base_index,6);
alpha_i = parameter_set(base_index,7);
zeta = epsilon*n0/kD;
sigma = R^2/d^2;
R_vals = parameter_set(slice,4)';
d_vals = parameter_set(slice,3)';
kappa_vals = parameter_set(slice, 6)';
zeta_vals = parameter_set(slice,1).*parameter_set(slice,2)./parameter_set(slice,5);

changing_value = 3;
param_1_vals = R^2./parameter_set(slice,changing_value).^2;
% param_1_vals = parameter_set(slice,changing_value);
E_all = E_all_total(:, slice);
alpha_A_vals = alpha_A_vals_total(slice);
alpha_B_vals = alpha_B_vals_total(slice);
phi_vals = phi_vals_total(slice);
h_phi_vals = h_phi_vals_total(slice);
S_A_vals = S_A_vals_total(slice);
S_B_vals = S_B_vals_total(slice);
Sigma_vals = Sigma_vals_total(slice);
sigma_vals = R_vals.^2./d_vals.^2;
lambda_vals = sqrt(kappa_vals./(kD*alpha_B_vals));

solution_set = solution_set_total(slice);

values_structure.E_all = E_all_total(:, slice);
values_structure.alpha_A_vals = alpha_A_vals_total(slice);
values_structure.alpha_B_vals = alpha_B_vals_total(slice);
values_structure.phi_vals = phi_vals_total(slice);
values_structure.h_phi_vals = h_phi_vals_total(slice);
values_structure.S_A_vals = S_A_vals_total(slice);
values_structure.S_B_vals = S_B_vals_total(slice);
values_structure.Sigma_vals = Sigma_vals_total(slice);
values_structure.sigma_vals = R_vals.^2./d_vals.^2;
values_structure.lambda_vals = sqrt(kappa_vals./(kD*alpha_B_vals));

values_structure.R_vals = parameter_set(slice,4)';
values_structure.d_vals = parameter_set(slice,3)';
values_structure.kappa_vals = parameter_set(slice, 6)';
values_structure.zeta_vals = parameter_set(slice,1).*parameter_set(slice,2)./parameter_set(slice,5);
values_structure.param_1_vals = param_1_vals;

values_structure.epsilon = parameter_set(base_index,1);
values_structure.n0 = parameter_set(base_index,2);
values_structure.d = parameter_set(base_index,3);
values_structure.R = parameter_set(base_index,4);
values_structure.kD = parameter_set(base_index,5);
values_structure.kappa = parameter_set(base_index,6);
values_structure.alpha_i = parameter_set(base_index,7);
values_structure.zeta = epsilon*n0/kD;
values_structure.sigma = R^2/d^2;

colours_cb = ["#27E0D8", "#648FFF", "#1D0C6B", "#DC267F", "#FE6100", ...
    "#FFB000", "#016D24"];

% close all

% xaxis_name = '$\alpha_i$';
xaxis_name = '$\sigma$';
% xaxis_name = '$R$ ($\mu$m)';
% xaxis_name = '$\kappa$ (pJ)';

% xscale = 'linear';
xscale = 'log';

anno_string = [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
    sprintf('$\\kappa = %0.2g$ pJ \n', kappa),...
    sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
    sprintf('$R = %0.2g$ $\\mu$m \n', R),...
    sprintf('$\\zeta = %0.2g$ \n', zeta)];

% patch_lims =  [1e-3, 1e-1];
patch_lims =  0;
% specific curves to plot
% plot_curves = [1,5, 10, 15, 25, 30, 35];
% plot_curves = round(linspace(15,70,7));
% plot_curves = round(linspace(1,24,5));
plot_curves = [20];
% plot_curves = [1,54,74,76,82,85,88];

% %% plotting absolute value

% E_absolute_value(xaxis_name, xscale, '$|E|$ (pJ)', 'log', param_1_vals, E_all,...
%     patch_lims, plot_curves, anno_string)
% 
% Energy_Subplots(xaxis_name, xscale, '$E$ (pJ)', 'lin', values_structure,...
%     patch_lims, plot_curves, anno_string)
% 
% E_value(xaxis_name, xscale, '$E$ (pJ)', 'lin', param_1_vals, E_all,...
%     patch_lims, plot_curves, anno_string)
% 
% dimensionless_parameters(xaxis_name, 'log', 'params', 'log',...
%     values_structure)
% 
% Non_Dim_Energy_All(xaxis_name, xscale, '$E/d^2 k_D$', 'lin', values_structure,...
%     patch_lims, plot_curves, anno_string)
% 
% Non_Dim_Energy_Separate(xaxis_name, xscale, '$E$ (pJ)', 'lin', values_structure,...
%     patch_lims, plot_curves, anno_string)
% 
% Non_Dim_Energy_Differences(xaxis_name, xscale, 'useless', 'lin', values_structure,...
%     patch_lims, plot_curves, anno_string)
% 
% Non_Dim_Separated_Energies(xaxis_name, xscale, '$E$', 'log',...
%     values_structure, patch_lims, plot_curves, 6)

% positive_percentages(xaxis_name, xscale, 'fraction of energy', 'lin',...
%     values_structure, patch_lims, plot_curves, 6)

% areas(xaxis_name, xscale, '$\%$ positive energy', 'lin',...
%     values_structure, patch_lims, plot_curves)
% 
% total_lipids(xaxis_name, xscale, '$\%$ positive energy', 'lin',...
%     values_structure, patch_lims, plot_curves)

if exist('fig_stretch', 'var') && isvalid(fig_stretch)
    fig_stretch = stretches(xaxis_name, xscale, '$\alpha_B$', 'lin',...
    values_structure, 0, 0, 0, fig_stretch);
else
    fig_stretch = stretches(xaxis_name, xscale, '$\alpha_B$', 'lin',...
        values_structure, 0, 0, 0);
end

if exist('fig_phi', 'var') && isvalid(fig_phi)
    fig_phi = phi_and_h_phi(xaxis_name, xscale, '$\phi$', 'lin',...
    values_structure, 0, 0,fig_phi);
else
    fig_phi = phi_and_h_phi(xaxis_name, xscale, '$\phi$', 'lin',...
    values_structure, 0, 0);
end

% energy_per_lipid(xaxis_name, xscale, '$E$', 'lin',...
%     values_structure, patch_lims, plot_curves)

% membrane_shapes(xaxis_name, xscale, '$E$', 'lin',...
%     values_structure, patch_lims, plot_curves)

% solution_space(xaxis_name, xscale, '$E$', 'lin',...
%     values_structure, patch_lims, solution_set, 11)

%%
% newcolors = ["#FFB000","#27E0D8","#1D0C6B","#DC267F","#FE6100","#648FFF","#016D24"];
% newcolors = ["#016D24","#648FFF","#FE6100","#DC267F","#1D0C6B","#27E0D8","#FFB000"];
% newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];
newcolors = flip(["#016D24","#FE6100","#DC267F"]);
% figure();
% plot(1:7,ones(1,7)'*[1:7])
colororder(newcolors)
