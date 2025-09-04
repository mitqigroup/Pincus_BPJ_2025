function fig_handle = plot_phi_minima(slices, data_structure, phi_vals, varargin)
    
% parse inputs and set defaults, starting at 4
args = varargin;
nargs = numel(args);
x_limits = nan;
% y_limits = nan;
use_lims = false;
plot_local_minima = false;
useLocalMinName = false;
plot_stars = false;
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
%     elseif strcmpi(args{k},'Ylimit')||strcmpi(args{k},'Y_limit')
%         y_limits = args{k+1};
%         use_lims = true;
%         k = k+1;
    elseif strcmpi(args{k},'PlotLocalMinima')||strcmpi(args{k},'Plot_Local_Minima')
        plot_local_minima = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'PlotStars')||strcmpi(args{k},'Plot_Stars')
        plot_stars = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'ColourOrder')||strcmpi(args{k},'ColorOrder')
        colours = args{k+1};
        % should be an array of strings as above
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
    case 1
        anno_string = [sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
            sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
            sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
            sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
        xlabel('$-\kappa/\epsilon n_0$')
    case 2
        anno_string = [sprintf('$\\epsilon n_0 = %0.2g$  pJ/$\\mu$m$^2$ \n', ps.epsilon_slice(1)),...
            sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
            sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
            sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
        xlabel('$-\kappa/\epsilon n_0$')
    case 3
        anno_string = [sprintf('$k_D = %0.2g$  pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
            sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
            sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
            sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
        xlabel('$\kappa/k_D$')
    case 5
        anno_string = [sprintf('$\\epsilon n_0 = %0.2g$  pJ/$\\mu$m$^2$ \n', ps.epsilon_slice(1)),...
            sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
            sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
            sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
        xlabel('$-k_D/\epsilon n_0$')
    case 6
        anno_string = [sprintf('$k_D = %0.2g$  pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
            sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
            sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
            sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
        xlabel('$-k_D/\epsilon n_0$')
    case 7
        anno_string = [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
            sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
            sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
            sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
            sprintf('$\\zeta = %0.2g$ \n', ps.zeta_slice(1))];
        xlabel('$\sigma$')
        useLocalMinName = true;
        
end

annotation('textbox', [0.185,0.671,0.257,0.207], 'String', anno_string)
ylabel('$\phi$')
axes1.XScale = 'log';
yticks([0, pi/4, pi/2, 3*pi/4, pi])
yticklabels({'$0$', '$\pi/4$', '$\pi/2$', '$3 \pi/4$', '$\pi$'})
ylim([-0.05,pi+0.1])
% legend('location', 'west')
legend('location', 'best')
colororder(newcolors)

% [~,~] = plot_cases_phi(ps, outputs_analyse, type);

%     % other slices
for ii=1:length(slices)
    if ii~=1
        slice = slices{ii};
        outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
            'PlotCurves', false, 'IgnoreNearNaN', false, 'CleanData', true);
        ps = get_param_slices(ds,slice);
    end
    if plot_stars
        [~,~] = plot_cases_phi(ps, outputs_analyse, type, 'Color', newcolors(ii), ...
            'PlotLocalMinima',plot_local_minima, 'UseLocalMinName', useLocalMinName,...
            'PlotStarsLimit', max(phi_vals)*0.999);
    else
        [~,~] = plot_cases_phi(ps, outputs_analyse, type, 'Color', newcolors(ii), ...
            'PlotLocalMinima',plot_local_minima, 'UseLocalMinName', useLocalMinName);
    end


end

if use_lims
    xlim(x_limits)
end

%     slice = slices{3};
%     outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
%         'PlotCurves', false, 'IgnoreNearNaN', false);
%     ps = get_param_slices(ds,slice);
%     [~,~] = plot_cases_phi(ps, outputs_analyse, type);


end