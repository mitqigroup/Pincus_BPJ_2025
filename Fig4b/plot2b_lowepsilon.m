clear variables

current_location = './data_epsilonpt3'; %specify file
ds = read_and_unpermute(current_location, 48, 'SortValues', true);

%%

% phi_vals = flip(deg2rad(linspace(0.01,179.9,150)));
% slice = [1,10,21,30,40,44,46];
slice = [1,11,30,39,40,45,47];
% slice = 39;

% newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];
newcolors = ["#003f5c", "#374c80", "#7a5195", "#bc5090", "#ef5675", "#ff764a", "#ffa600"];

ps = get_param_slices(ds, slice);

% %% analyse data from each, we want local minima and maxima and end points.
E_input = ds.E_all_nonlinear;
% E_input(ds.E_all_toroid<ds.E_all_nonlinear) = nan;
outputs_analyse = analyse_phi_curve(E_input,ds, slice, ds.phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false, 'RemoveUnphysical', true);
% outputs_analyse = analyse_phi_curve(E_input,ds, slice, phi_vals,...
%     'PlotCurves', true, 'IgnoreNearNaN', false, 'RemoveUnphysical', true);
% outputs_analyse = analyse_phi_curve(E_input,ds, 1020:1030, phi_vals,...
%     'PlotCurves', true, 'IgnoreNearNaN', false);


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
% f1 = plot_phi_curves(ds, slice, phi_vals,...
%     'AnalysisData', outputs_analyse,'PlotMinima', true,...
%     'RelativeEnergy', true,'ylabel','$\Delta E/\mathrm{max}(|\Delta E|)$',...
%     'DisplayNameIndex',8,'AnnotationString', anno_string,...
%     'ColourOrder', newcolors, 'ThickerLines', []);
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

% %% Sigma vs d
% slice = 1:48;
% Sig_max = ds.Sigma_vals_nonlinear(end,:);
% ps = get_param_slices(ds, slice);
% E_input = ds.E_all_nonlinear;
% outputs_analyse = analyse_phi_curve(E_input,ds, slice, ds.phi_vals,...
%     'PlotCurves', false, 'IgnoreNearNaN', false, 'RemoveUnphysical', true);
% Sig_i = (ps.alpha_i_slice.*ps.kD_slice)';
% figure();
% hold on
% axes1 = gca;
% axes1.XScale = 'log';
% axes1.YScale = 'log';
% plot(ps.d_slice./ps.R_slice, outputs_analyse.Sigma_at_min./Sig_i)
% plot(ps.d_slice./ps.R_slice, Sig_max./Sig_i)
% plot(ps.d_slice./ps.R_slice, (Sig_i+4*ps.sigma_slice'.*ps.kD_slice')./Sig_i)
% plot(ps.d_slice./ps.R_slice, ones(size(ps.d_slice)), ':')
% % plot(ps.d_slice./ps.R_slice, ps.alpha_i_slice.*ps.kD_slice, ':')
% 
% %% construct slices, can't really be done in a separate function
% % generate list of independent variables to run, which should be in the
% % order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
% R_vals = [0.05];                 % um
% sigma_vals = [logspace(-9,log10(2e-6),10),linspace(3e-6,4.5e-6,28),logspace(log10(5e-6),-1,10)];
% %d_vals = sqrt(R_vals.^2./sigma_vals);           % um
% kD_vals = 300/10^12*1e9;                        % picoJ/um^2
% kD_base = 300/10^12*1e9;                        % picoJ/um^2
% zeta_vals = 0.0002848;                  % dimensionless
% epsilon_vals = -zeta_vals*kD_base;              % picoJ/um^2
% n0_vals = 1;                                    % fraction
% kappa_vals = 1e-19*1e12;         % picoJ
% % kappa_vals = logspace(-21,-15, 60)*1e12;
% alpha_i_vals = [1e-5];                            % fraction
% phi_vals = flip(deg2rad(linspace(0.01,179.9,150)));
% % other constants
% color_number = 5;
% 
% ii = 0;
% slice = [];
% for rr = 1:length(R_vals)
%     for ss = 1:length(sigma_vals)
%         for kk = 1:length(kD_vals)
%             for ee = 1:length(epsilon_vals)
%                 for nn = 1:length(n0_vals)
%                     for pp = 1:length(kappa_vals)
%                         for aa = 1:length(alpha_i_vals)
%                             ii = ii+1;
% %                             if any(ss==[1,10,21,30,40,44,46])
% %                                 slice = [slice, ii];
% %                             end
%                             if true
%                                 slice = [slice, ii];
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% size_params = size(ds.parameter_set);
% base_index = slice(1);
