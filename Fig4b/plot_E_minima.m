function fig_handle = plot_E_minima(slices, data_structure, phi_vals, varargin)


% parse inputs and set defaults, starting at 4
args = varargin;
nargs = numel(args);
x_limits = nan;
y_limits = nan;
use_lims = false;
plot_specific_points = false;
type = 1;
k = 1;
while k<=nargs
    % if we get a figure handle, use it
    if strcmpi(args{k},'type')
        type = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'Xlimit')||strcmpi(args{k},'X_limit')
        x_limits = args{k+1};
        use_lims = true;
        k = k+1;
    elseif strcmpi(args{k},'Ylimit')||strcmpi(args{k},'Y_limit')
        y_limits = args{k+1};
        use_lims = true;
        k = k+1;
    elseif strcmpi(args{k},'PlotSpecificPoints')||strcmpi(args{k},'Plot_Specific_Points')
        plotting_points = args{k+1};
        plot_specific_points = true;
        % this should be a cell array of specific points to plot, with the
        % same size as the slices
        assert(length(slices)==length(plotting_points), 'Slices not same as plotting points');
        k = k+1;
    else
        error('intput not recognised!');
    end
    k = k+1;
end
    
newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];

ds = data_structure;

slice = slices{1};
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false, 'RemoveUnphysical', true);

ps = get_param_slices(ds,slice);

fig_handle = figure('Position',[400,100,800,600]);
axes1 = gca;
hold on

switch type
    case 4
        anno_string = [sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
            sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
            sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
        sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
    xlabel('$\kappa/k_D$')
    ylabel('$E_\mathrm{bend}/E_\mathrm{stretch}$')
case 1
    anno_string = [sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
        sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
        sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
        sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
    xlabel('$-\kappa/\epsilon n_0$')
    ylabel('-$E_\mathrm{bend}/E_\mathrm{adhesion}$')
case 2
    anno_string = [sprintf('$\\epsilon n_0 = %0.2g$  pJ/$\\mu$m$^2$  \n', ps.epsilon_slice(1)),...
        sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
        sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
        sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
    xlabel('$-\kappa/\epsilon n_0$')
    ylabel('-$E_\mathrm{bend}/E_\mathrm{adhesion}$')
case 3
    anno_string = [sprintf('$k_D = %0.2g$  pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
        sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
        sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
        sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
    xlabel('$\kappa/k_D$')
    ylabel('$E_\mathrm{bend}/\Delta E_\mathrm{stretch}$')
case 5
    anno_string = [sprintf('$\\epsilon n_0 = %0.2g$  pJ/$\\mu$m$^2$ \n', ps.epsilon_slice(1)),...
        sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
        sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
        sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
    xlabel('$-k_D/\epsilon n_0$')
    ylabel('$-\Delta E_\mathrm{stretch}/E_\mathrm{adhesion}$')
case 6
    anno_string = [sprintf('$k_D = %0.2g$  pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
        sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
        sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
        sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
    xlabel('$-k_D/\epsilon n_0$')
    ylabel('$-\Delta E_\mathrm{stretch}/E_\mathrm{adhesion}$')

end

ps = get_param_slices(ds,slice);
annotation('textbox', [0.646,0.152,0.247,0.21], 'String', anno_string);
axes1.XScale = 'log';
axes1.YScale = 'log';
if use_lims
    ylim(y_limits)
    xlim(x_limits)
end
[~,~] = plot_cases_E(ps, outputs_analyse, type);
legend('location', 'best')
colororder(newcolors)

for ii=2:length(slices)
    slice = slices{ii};
    outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
        'PlotCurves', false, 'IgnoreNearNaN', false);
    ps = get_param_slices(ds,slice);
    [~,~] = plot_cases_E(ps, outputs_analyse, type);
end
    
%     slice = slices{3};
%     outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
%         'PlotCurves', false, 'IgnoreNearNaN', false);
%     ps = get_param_slices(ds,slice);
%     [~,~] = plot_cases_E(ps, outputs_analyse, type);

end