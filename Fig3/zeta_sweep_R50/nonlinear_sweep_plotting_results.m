clear variables

current_location = ['D:\Work\Supercloud\membranes\',...
    'nonlinear_parameter_sweeps\ep_sig_sweep\zeta_sweep_R50\middle'];
ds_middle = read_and_unpermute(current_location, 48, 'sortvalues', true);
current_location = ['D:\Work\Supercloud\membranes\',...
    'nonlinear_parameter_sweeps\ep_sig_sweep\zeta_sweep_R50\down'];
ds_down = read_and_unpermute(current_location, 48, 'sortvalues', true);
newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];


%% smaller tests
% rv = 1;
% sv = 1;
% kv = 1;
% pv = 1;
% av = 1;
% ev = 22;
% slice = get_slice(rv,sv,kv,ev, pv, av);
% outputs_analyse = analyse_phi_curve_updown(ds_middle,ds_down,slice, 'plotcurves', true);
close all
slice = [1, 49, 130, 192,48];
ps = get_param_slices(ds_middle,slice);

% %% plot against phi
anno_string = [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
    sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
    sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
    sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
    sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
f1 = plot_phi_curves_updown(ds_middle, ds_down, slice,...
    'PlotMinima', true,'RemoveUnphysical', true,...
    'RelativeEnergy', true,'ylabel','$\Delta E/\mathrm{max}(|\Delta E|)$',...
    'DisplayNameIndex',9,'AnnotationString', anno_string,...
    'ColourOrder', newcolors, 'ThickerLines', [1,5], 'PlotLocalMinima', true);
f1.Children(1).Location='northwest';
plot([0,180],[-0.213821,-0.213821],':', 'Color', newcolors(3), 'HandleVisibility','off')

% %%
slice_new = 1:192;
type = 8;
f1 = plot_phi_minima_updown(ds_middle, ds_down,{slice_new}, 'type', type,...
    'Plot_Local_Minima', true, 'PlotStars', true,'PlotVerticalLines', slice);
% change colours of main lines to black?
f1.Children(2).Children(1).Color = 'k';
f1.Children(2).Children(2).Color = 'k';
f1.Children(2).Children(1).LineWidth = 2;
f1.Children(2).Children(2).LineWidth = 2;
% xlim([2.5e-4, 2.9e-4])
xlim([5e-6, 2e-3])

%% zoomed in
% rv = 1;
% sv = 1;
% kv = 1;
% pv = 1;
% av = 1;
% ev = 22;
% slice = get_slice(rv,sv,kv,ev, pv, av);
% outputs_analyse = analyse_phi_curve_updown(ds_middle,ds_down,slice, 'plotcurves', true);
% close all
slice = [49, 60, 100,130, 160];
ps = get_param_slices(ds_middle,slice);
% outputs_analyse = analyse_phi_curve_updown(ds_middle,ds_down,slice, 'plotcurves', true);

% plot against phi
anno_string = [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
    sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
    sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
    sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
    sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];
f1 = plot_phi_curves_updown(ds_middle, ds_down, slice,...
    'PlotMinima', true,'RemoveUnphysical', true,...
    'RelativeEnergy', true,'ylabel','$\Delta E/\mathrm{max}(|\Delta E|)$',...
    'DisplayNameIndex',9,'AnnotationString', anno_string,...
    'ColourOrder', newcolors, 'ThickerLines', [], 'PlotLocalMinima', true);
f1.Children(1).Location='northwest';
plot([0,180],[-0.213821,-0.213821],':', 'Color', newcolors(4), 'HandleVisibility','off')

% %%
slice_new = 1:192;
type = 8;
f1 = plot_phi_minima_updown(ds_middle, ds_down,{slice_new}, 'type', type,...
    'Plot_Local_Minima', true, 'PlotStars', true,'PlotVerticalLines', slice);
% change colours of main lines to black?
f1.Children(2).Children(1).Color = 'k';
f1.Children(2).Children(2).Color = 'k';
f1.Children(2).Children(1).LineWidth = 2;
f1.Children(2).Children(2).LineWidth = 2;
xlim([2.63e-4, 2.8e-4])
% xlim([5e-6, 2e-3])


%% functions

function slice = get_slice(rv, sv, kv, ev, pv, av)
R_vals = [0.05];                 % um
sigma_vals = 3.1e-06;                % surface fraction
kD_vals = 300/10^12*1e9;                        % picoJ/um^2
kD_base = 300/10^12*1e9;                        % picoJ/um^2
zeta_vals = [logspace(-5,-3,48),linspace(0.0002653, 0.00028169, 48*3)];                  % dimensionless
epsilon_vals = -zeta_vals*kD_base;              % picoJ/um^2
n0_vals = 1;                                    % fraction
kappa_vals = 1e-19*1e12;         % picoJ
alpha_i_vals = [1e-6];                            % fraction

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
                            if any(rr==rv) && any(ss==sv) && any(kk==kv) && any(ee==ev) && any(rr==rv) && any(pp==pv) && any(aa==av)
                                slice = [slice, ii];
                            end
                        end
                    end
                end
            end
        end
    end
end

end











