clear variables

current_location = ['Z:\Isaac\Isaac_Data\Supercloud\membranes\',...
    'nonlinear_parameter_sweeps\ep_sig_sweep\example_1_sig_sweep\'];
% main data structure
ds = read_and_unpermute(current_location, 48);

%% analyse data from each, we want local minima and maxima and end points.
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals, 0);

%% plot against phi
f1 = plot_phi_curves(ds, slice, phi_vals,...
    'AnalysisData', outputs_analyse,'PlotMinima', true, 'ColorOrder', newcolors);

%% kappa constant, changing epsilon, three k_D values, phi plot

phi_vals = flip(deg2rad(linspace(0.01,179.9,150)));

newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];

% kD = 0.01
slice = 8:24:560;
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false);

ps = get_param_slices(ds,slice);

anno_string = [sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
    sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
    sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
    sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];

f3 = figure('Position',[400,100,800,600]);
annotation('textbox', 'String', anno_string)
axes1 = gca;
hold on
xlabel('$-\kappa/\epsilon n_0$')
ylabel('$\phi$')
axes1.XScale = 'log';
yticks([0, pi/4, pi/2, 3*pi/4, pi])
yticklabels({'$0$', '$\pi/4$', '$\pi/2$', '$3 \pi/4$', '$\pi$'})
xticks(logspace(-5,-2,4))
ylim([-0.05,pi+0.1])
xlim([2e-5, 1e-2])
x = -ps.kappa_slice./(ps.epsilon_slice.*ps.n0_slice);
y = outputs_analyse.phi_at_min;
v1 = plot(x, y,'displayname', sprintf('$k_D = %0.3g$', ps.kD_slice(1)));
legend('location', 'sw')
colororder(newcolors)

% kD = 0.1
slice = 584:24:1136;
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false);

ps = get_param_slices(ds,slice);
x = -ps.kappa_slice./(ps.epsilon_slice.*ps.n0_slice);
y = outputs_analyse.phi_at_min;
v1 = plot(x, y,'displayname', sprintf('$k_D = %0.3g$', ps.kD_slice(1)));

% kD = 1
slice = 1160:24:1712;
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false);

ps = get_param_slices(ds,slice);
x = -ps.kappa_slice./(ps.epsilon_slice.*ps.n0_slice);
y = outputs_analyse.phi_at_min;
v1 = plot(x, y,'displayname', sprintf('$k_D = %0.3g$', ps.kD_slice(1)));

%% kappa constant, changing epsilon, three k_D values, Ebend/E_adhesion plot

phi_vals = flip(deg2rad(linspace(0.01,179.9,150)));

newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];

% kD = 0.01
slice = 8:24:560;
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false);

ps = get_param_slices(ds,slice);
E_unwrapped = ps.kD_slice./2.*ps.alpha_i_slice.^2.*ps.d_slice.^2./(1+ps.alpha_i_slice);

anno_string = [sprintf('$\\kappa = %0.2g$ pJ \n', ps.kappa_slice(1)),...
    sprintf('$\\alpha_i = %0.2g$ \n', ps.alpha_i_slice(1)),...
    sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
    sprintf('$\\sigma = %0.2g$ \n', ps.sigma_slice(1))];

f3 = figure('Position',[400,100,800,600]);
annotation('textbox', 'String', anno_string)
axes1 = gca;
hold on
xlabel('$-\kappa/\epsilon n_0$')
ylabel('$-E_\mathrm{bend}/E_\mathrm{adhesion}$')
axes1.XScale = 'log';
axes1.YScale = 'log';
xticks(logspace(-5,-2,4))
xlim([2e-5, 1e-2])
x = -ps.kappa_slice./(ps.epsilon_slice.*ps.n0_slice);
y = -(sum(outputs_analyse.E_all_min(5:6,:),1))./outputs_analyse.E_all_min(2,:);
y(outputs_analyse.phi_at_min<0.01) = nan;
v1 = plot(x, y,'displayname', sprintf('$k_D = %0.3g$', ps.kD_slice(1)));
legend('location', 'sw')
colororder(newcolors)

% kD = 0.1
slice = 584:24:1136;
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false);

ps = get_param_slices(ds,slice);
x = -ps.kappa_slice./(ps.epsilon_slice.*ps.n0_slice);
y = -(sum(outputs_analyse.E_all_min(5:6,:),1))./outputs_analyse.E_all_min(2,:);
y(outputs_analyse.phi_at_min<0.01) = nan;
v1 = plot(x, y,'displayname', sprintf('$k_D = %0.3g$', ps.kD_slice(1)));

slice = 1160:24:1712;
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false);

ps = get_param_slices(ds,slice);

% kD = 1
x = -ps.kappa_slice./(ps.epsilon_slice.*ps.n0_slice);
y = -(sum(outputs_analyse.E_all_min(5:6,:),1))./outputs_analyse.E_all_min(2,:);
y(outputs_analyse.phi_at_min<0.01) = nan;
v1 = plot(x, y,'displayname', sprintf('$k_D = %0.3g$', ps.kD_slice(1)));


%% construct slices, can't really be done in a separate function
% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = [0.05];                 % um
sigma_vals = logspace(-9,-2,48);                % surface fraction
%d_vals = sqrt(R_vals.^2./sigma_vals);           % um
kD_vals = 300/10^12*1e9;                        % picoJ/um^2
kD_base = 300/10^12*1e9;                        % picoJ/um^2
zeta_vals = 0.0002848;                  % dimensionless
epsilon_vals = -zeta_vals*kD_base;              % picoJ/um^2
n0_vals = 1;                                    % fraction
kappa_vals = 1e-19*1e12;         % picoJ
% kappa_vals = logspace(-21,-15, 60)*1e12;
alpha_i_vals = [1e-5];                            % fraction
phi_vals = flip(deg2rad(linspace(0.01,179.9,150)));
% other constants
N = 3e5;
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
%                             if ss>1 && ss<30 && mod(ss,4)==0 && rr==1
%                                 slice = [slice, ii];
%                             end
%                             if (ss==1 || ss==15 || ss==21 || ss==22|| ss==23 || ss==24 || ss==26) && rr==1
%                                 slice = [slice, ii];
%                             end
                            if (ss==1 || ss==10 || ss==30 || ss==40 || ss==46) && rr==1
                                slice = [slice, ii];
                            end
                        end
                    end
                end
            end
        end
    end
end
size_params = size(ds.parameter_set);
base_index = slice(1);
newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];

% updated parameter slices
epsilon_slice = ds.parameter_set(slice,1);
n0_slice = ds.parameter_set(slice,2);
d_slice = ds.parameter_set(slice,3);
R_slice = ds.parameter_set(slice,4);
kD_slice = ds.parameter_set(slice,5);
kappa_slice = ds.parameter_set(slice,6);
alpha_i_slice = ds.parameter_set(slice,7);
sigma_slice = pi*R_slice.^2./d_slice.^2;
zeta_slice = epsilon_slice.*n0_slice./kD_slice;