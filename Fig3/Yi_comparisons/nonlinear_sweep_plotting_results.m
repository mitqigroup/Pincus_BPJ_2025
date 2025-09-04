clear variables
addpath("Finer_grained")

current_location = ['D:\Work\Supercloud\membranes\nonlinear_parameter_sweeps\Yi_comparisons\'];
% main data structure
ds = read_and_unpermute(current_location, 12);

% %% plot the Yi wrapping results here
%% plot against phi
newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];
slice = 1:6;
phi_vals = flip(deg2rad(linspace(0.01,179.9,150)));
ps = get_param_slices(ds,slice);
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false);
f1 = plot_phi_curves(ds, slice, phi_vals,...
    'ylabel','$\Delta E / \pi \kappa$',...
    'AnalysisData', outputs_analyse,'PlotMinima', true,...
    'PlotLocalMinima', true,'ColorOrder', newcolors, ...
    'RelativeEnergy', false, 'ScaleFactor', pi*ps.kappa_slice','DisplayNameIndex',10);
anno_string = [sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
    sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
    sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
    sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
    sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
annotation(f1,'textbox', [0.35,0.66,0.26,0.241],'String', anno_string)

% annotation(f1,'textbox', [0.35,0.66,0.26,0.241],'String',...
%     'global minima')
% annotation(f1,'textbox', [0.35,0.66,0.26,0.241],'String',...
%     'local minima')

%% plot against phi with Yi data
newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];


data = readmatrix('yi_data.csv');

size_data = size(data);

% figure();
hold on

for ii=1:2:size_data(2)
    plot(rad2deg(acos(1-2*data(:,ii))), data(:,ii+1), 'o',...
        'Color',newcolors((ii+1)/2), 'HandleVisibility','off');
end

annotation('textbox', 'string', [sprintf('$\\circ$ Yi \n'), sprintf('$-$ Our Work')]);

%% plot relative
newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];
slice = 1:6;
phi_vals = flip(deg2rad(linspace(0.01,179.9,150)));
ps = get_param_slices(ds,slice);
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false, 'RemoveUnphysical', true);
f1 = plot_phi_curves(ds, slice, phi_vals,...
    'ylabel','$\Delta E / \max(|\Delta E|)$',...
    'AnalysisData', outputs_analyse,'PlotMinima', true, 'RemoveUnphysical', true,...
    'PlotLocalMinima', true,'ColorOrder', newcolors, ...
    'RelativeEnergy', true, 'ScaleFactor', pi*ps.kappa_slice','DisplayNameIndex',10);
anno_string = [sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
    sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
    sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
    sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
    sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
annotation(f1,'textbox', [0.35,0.66,0.26,0.241],'String', anno_string)

% annotation(f1,'textbox', [0.35,0.66,0.26,0.241],'String',...
%     'global minima')
% annotation(f1,'textbox', [0.35,0.66,0.26,0.241],'String',...
%     'local minima')

%% plot against phi with Yi data relative
newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];


data = readmatrix('yi_data.csv');

size_data = size(data);

% figure();
hold on

for ii=1:2:size_data(2)
    plot(rad2deg(acos(1-2*data(:,ii))), data(:,ii+1)/max(abs(data(:,ii+1))), 'o',...
        'Color',newcolors((ii+1)/2), 'HandleVisibility','off');
end

annotation('textbox', 'string', [sprintf('$\\circ$ Yi \n'), sprintf('$-$ Our Work')]);

%% construct slices, can't really be done in a separate function
% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = [0.1];                 % um
sigma_vals = [1e-9,1e-5,1e-2];                % surface fraction
%d_vals = sqrt(R_vals.^2./sigma_vals);           % um
kD_vals = 0.3;                        % picoJ/um^2
% zeta_vals = 0.0002848;                  % dimensionless
epsilon_vals = -[2.2,3,4,5,6.5,7.7]*1e-7/(2*R_vals^2);              % picoJ/um^2
n0_vals = 1;                                    % fraction
kappa_vals = 1e-7;         % picoJ
% kappa_vals = logspace(-21,-15, 60)*1e12;
alpha_i_vals = 1e-7/0.3*2/(2*R_vals^2);                            % fraction
phi_vals = flip(deg2rad(linspace(0.01,179.9,150)));
% other constants
color_number = 5;

ii = 0;
slice = [];
for rr = 1:length(R_vals)
    for ss = 1:length(sigma_vals)
        for kk = 1:length(kD_vals)
            for ee = 1:length(epsilon_vals)
                for nn = 1:length(n0_vals)
                    for pp = 1:length(kappa_vals)
                        for aa = 1:length(alpha_i_vals)
                            ii = ii+1;
                            if ss==1
                                slice = [slice, ii];
                            end
                        end
                    end
                end
            end
        end
    end
end
% size_params = size(parameter_set);
% base_index = slice(1);
slice
