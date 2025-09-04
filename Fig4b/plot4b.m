clear variables
close all
%fig_handle = figure('Position',[400,100,800,600]);

newcolors = ["#003f5c", "#374c80", "#7a5195", "#bc5090", "#ef5675", "#ff764a", "#ffa600"];
type = 7;
epval=zeros(1,3);
hold on;

current_location = './data_epsilonpt3'; %specify file
ds = read_and_unpermute(current_location, 48, 'SortValues', true);
slice = 1:48;
ps = get_param_slices(ds, slice);
epval(1)=ps.epsilon_slice(1);
%% plot against phi
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, ds.phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false, 'RemoveUnphysical', true);
plot_cases_phi(ps, outputs_analyse, type, 'Color', newcolors(1), ...
            'PlotLocalMinima',true, 'UseLocalMinName', false,...
            'PlotStarsLimit', max(ds.phi_vals)*0.999);

current_location = './data_epsilon1'; %specify file
ds = read_and_unpermute(current_location, 1, 'SortValues', true);
slice=1:18;
ps = get_param_slices(ds, slice);
epval(2)=ps.epsilon_slice(1);

%% plot against phi
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, ds.phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false, 'RemoveUnphysical', true);
plot_cases_phi(ps, outputs_analyse, type, 'Color', newcolors(5), ...
            'PlotLocalMinima',true, 'UseLocalMinName', false,...
            'PlotStarsLimit', max(ds.phi_vals)*0.999);

current_location = './data_epsilon3'; %specify file
ds = read_and_unpermute(current_location, 1, 'SortValues', true);
slice=1:18;
ps = get_param_slices(ds, slice);
epval(3)=ps.epsilon_slice(1);

%% plot against phi
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, ds.phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false, 'RemoveUnphysical', true);

plot_cases_phi(ps, outputs_analyse, type, 'Color', newcolors(7), ...
            'PlotLocalMinima',true, 'UseLocalMinName', false,...
            'PlotStarsLimit', max(ds.phi_vals)*0.999);


hold off;
%annotation('textbox', [0.185,0.671,0.257,0.207], 'String', anno_string)
ylabel('$\phi$')
xlabel('$\sigma$')
axes1 = gca;
axes1.XScale = 'log';
yticks([0, pi/4, pi/2, 3*pi/4, pi])
yticklabels({'$0$', '$\pi/4$', '$\pi/2$', '$3 \pi/4$', '$\pi$'})
ylim([0.15*pi,0.8*pi])
xlim([0.1*1e-3,0.5*1e-1])
colororder(newcolors)
string=[];
for i=1:3;
string{i}=sprintf('$\\epsilon n_0= %0.2g$ mJ/m$^2$\n', epval(i)*1000);
end
legend(string{1},string{2},string{3});
anno_string = [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
    sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
    sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
    sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
    sprintf('$\\zeta = %0.2g$ \n', ps.zeta_slice(1)),...
    sprintf('$\\epsilon n_0= %0.2g$ mJ/m$^2$\n', ps.epsilon_slice(1)*1000)]


