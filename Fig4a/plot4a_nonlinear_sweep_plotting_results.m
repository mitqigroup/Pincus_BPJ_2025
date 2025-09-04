clear variables

current_location = '.'; %specify file
ds = read_and_unpermute(current_location, 48, 'SortValues', true);

slice = [1,11,30,39,40,45,47];

% newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];
newcolors = ["#003f5c", "#374c80", "#7a5195", "#bc5090", "#ef5675", "#ff764a", "#ffa600"];

ps = get_param_slices(ds, slice);

% %% analyse data from each, we want local minima and maxima and end points.
E_input = ds.E_all_nonlinear;
outputs_analyse = analyse_phi_curve(E_input,ds, slice, ds.phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false, 'RemoveUnphysical', true);

%%
counter = 0;
for index=slice
    counter = counter+1;
    R = ds.parameter_set(index,4);
    d = ds.parameter_set(index,3);
    phi = outputs_analyse.phi_at_min(counter);
    kappa = ds.parameter_set(index,6);
    Sigma = outputs_analyse.Sigma_at_min(counter);
    plot_shape_no_axis(R, d, phi, kappa, Sigma, 'colour', newcolors(counter))
end

%% plot against phi
anno_string = [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
    sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
    sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
    sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
    sprintf('$\\zeta = %0.2g$ \n', ps.zeta_slice(1))];
f1 = plot_phi_curves(ds, slice, ds.phi_vals,...
    'AnalysisData', outputs_analyse,'PlotMinima', true,...
    'RelativeEnergy', true,'ylabel','$\Delta E/\mathrm{max}(|\Delta E|)$',...
    'DisplayNameIndex',8,'AnnotationString', anno_string,...
    'ColourOrder', newcolors, 'ThickerLines', [],...
    'RemoveUnphysical', true,'PlotLocalMinima',true); %2(a) bigger plot
annotation('rectangle',[0.131,0.143,0.50325,0.415],LineStyle='--',LineWidth=0.5)
plot([0,180],[-1,-1],':', 'Color', newcolors(3), 'HandleVisibility','off')
f1.Children(1).String{1} = '$C, \sigma \rightarrow 0$';
f1.Children(1).Location='northwest';
f2 = plot_phi_curves(ds, slice, ds.phi_vals,...
    'AnalysisData', outputs_analyse,'PlotMinima', true,...
    'RelativeEnergy', false,'ylabel','$\Delta E$',...
    'DisplayNameIndex',8,'IncludeLegend', false,...
    'ColourOrder', newcolors, 'ThickerLines', [],'PlotLocalMinima',true);
plot([0,180],[-2.36904e-6,-2.36904e-6],':', 'Color', newcolors(3), 'HandleVisibility','off')
xlim([0,120]);
ylim([-5e-6,2e-6]); %inset

%% new
slices{1} = 1:48;
type = 7;
x_limits = [1e-9,1e-1];
f1 = plot_phi_minima(slices, ds, ds.phi_vals, 'type', type,...
    'Plot_Local_Minima', true,'PlotStars',true, 'XLimit', x_limits);%part (b)?

