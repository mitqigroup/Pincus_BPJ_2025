clear variables
close all

[filePath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
current_location = [filePath,'/middle'];
ds_middle = read_and_unpermute(current_location, 48, 'sortvalues', true);
current_location = [filePath,'/down'];
ds_down = read_and_unpermute(current_location, 48, 'sortvalues', true);

newcolors = ["#003f5c", "#374c80", "#7a5195", "#bc5090", "#ef5675", "#ff764a", "#ffa600"];


%% lambda plots

%sv = 1:192;
sv = 1:1000;
kv = 1;
ev = 1;
pv = 1;
rv = 1; %50 nm
av = 2;

%figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% axes1.YScale = 'log';


slice = get_slice(rv,sv,kv,ev, pv, av);
outputs_analyse = analyse_phi_curve_updown(ds_middle,ds_down,slice, ...
    'MovmedianRange', 5,'Movmedian_Threshold',0.1,'NSplinePoints',1000);
ps = get_param_slices(ds_middle,slice);
alpha_i = ps.alpha_i_slice(1);
kappa = ps.kappa_slice(1);
kD = ps.kD_slice(1);
R = ps.R_slice(1);
d = ps.d_slice';
S = [cell2mat(outputs_analyse.all_global_minima).S_A_min] + ...
    [cell2mat(outputs_analyse.all_global_minima).S_B_min];
alpha_B = [cell2mat(outputs_analyse.all_global_minima).alpha_B_min];
Sigma = kD*(S./d.^2*(1+alpha_i)-1);
lambda = sqrt(kappa./Sigma);
Sigma_B = kD*alpha_B;
lambda_B = sqrt(kappa./Sigma_B);
Sig_i = kD*alpha_i;
plot(d/R, log(Sigma/Sig_i), 'color', newcolors(1), ...
    'DisplayName',sprintf('$\\Sigma_i = %0.2g$ $\\mu$J/m$^2$, $R = %0.2g$ nm', alpha_i*kD*10^6,R*1000))
%plot(d/R, Sigma_B/Sig_i, ':', 'color', newcolors(1), 'HandleVisibility','off')

av=3;
slice = get_slice(rv,sv,kv,ev, pv, av);
outputs_analyse = analyse_phi_curve_updown(ds_middle,ds_down,slice, ...
    'MovmedianRange', 5,'Movmedian_Threshold',0.1,'NSplinePoints',1000);
ps = get_param_slices(ds_middle,slice);
alpha_i = ps.alpha_i_slice(1);
kappa = ps.kappa_slice(1);
kD = ps.kD_slice(1);
R = ps.R_slice(1);
d = ps.d_slice';
S = [cell2mat(outputs_analyse.all_global_minima).S_A_min] + ...
    [cell2mat(outputs_analyse.all_global_minima).S_B_min];
alpha_B = [cell2mat(outputs_analyse.all_global_minima).alpha_B_min];
Sigma = kD*(S./d.^2*(1+alpha_i)-1);
lambda = sqrt(kappa./Sigma);
Sigma_B = kD*alpha_B;
lambda_B = sqrt(kappa./Sigma_B);
Sig_i = kD*alpha_i;
plot(d/R, log(Sigma/Sig_i), 'color', newcolors(5), ...
    'DisplayName',sprintf('$\\Sigma_i = %0.2g$ $\\mu$J/m$^2$, $R = %0.2g$ nm', alpha_i*kD*10^6,R*1000))

% av=2;
% rv=2;
% slice = get_slice(rv,sv,kv,ev, pv, av);
% outputs_analyse = analyse_phi_curve_updown(ds_middle,ds_down,slice, ...
%     'MovmedianRange', 5,'Movmedian_Threshold',0.1,'NSplinePoints',1000);
% ps = get_param_slices(ds_middle,slice);
% alpha_i = ps.alpha_i_slice(1);
% kappa = ps.kappa_slice(1);
% kD = ps.kD_slice(1);
% R = ps.R_slice(1);
% d = ps.d_slice';
% S = [cell2mat(outputs_analyse.all_global_minima).S_A_min] + ...
%     [cell2mat(outputs_analyse.all_global_minima).S_B_min];
% alpha_B = [cell2mat(outputs_analyse.all_global_minima).alpha_B_min];
% Sigma = kD*(S./d.^2*(1+alpha_i)-1);
% lambda = sqrt(kappa./Sigma);
% Sigma_B = kD*alpha_B;
% lambda_B = sqrt(kappa./Sigma_B);
% Sig_i = kD*alpha_i;
% plot(d/R, log(Sigma/Sig_i), 'color', newcolors(7), ...
%     'DisplayName',sprintf('$\\Sigma_i = %0.2g$ $\\mu$J/m$^2$, $R = %0.3g$ nm', alpha_i*kD*10^6,R*1000))

av=3;
rv=2;
slice = get_slice(rv,sv,kv,ev, pv, av);
outputs_analyse = analyse_phi_curve_updown(ds_middle,ds_down,slice, ...
    'MovmedianRange', 5,'Movmedian_Threshold',0.1,'NSplinePoints',1000);
ps = get_param_slices(ds_middle,slice);
alpha_i = ps.alpha_i_slice(1);
kappa = ps.kappa_slice(1);
kD = ps.kD_slice(1);
R = ps.R_slice(1);
d = ps.d_slice';
S = [cell2mat(outputs_analyse.all_global_minima).S_A_min] + ...
    [cell2mat(outputs_analyse.all_global_minima).S_B_min];
alpha_B = [cell2mat(outputs_analyse.all_global_minima).alpha_B_min];
Sigma = kD*(S./d.^2*(1+alpha_i)-1);
lambda = sqrt(kappa./Sigma);
Sigma_B = kD*alpha_B;
lambda_B = sqrt(kappa./Sigma_B);
Sig_i = kD*alpha_i;
plot(d/R, log(Sigma/Sig_i), 'color', newcolors(7), ...
    'DisplayName',sprintf('$\\Sigma_i = %0.2g$ $\\mu$J/m$^2$, $R = %0.3g$ nm', alpha_i*kD*10^6,R*1000))

plot([3,100], [1.1 1.1], 'k-.','DisplayName','$d_c/R$')
ylim([1,7]);
xlim([3,100]);
xlabel('$d/R$')
ylabel('$\Sigma/\Sigma_i$')

legend

% a3 = annotation('textbox', 'String',...
%     [sprintf('$\\epsilon n_0 = %0.3g$ mJ/m$^2$ \n', ps.epsilon_slice(1)*1000),...
%     sprintf('$\\Sigma_i = %0.3g$ \n', ps.alpha_i_slice(1)*ps.kD_slice(1)),...
%     sprintf('$\\kappa = %0.3g$ J \n', ps.kappa_slice(1)),...
%      sprintf('$R = %0.2g$ $\\mu$m \n', ps.R_slice(1)),...
%     sprintf('$k_D = %0.3g$ N/m', ps.kD_slice(1))])

%% functions

function slice = get_slice(rv, sv, kv, ev, pv, av)
% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = [0.05,0.5];                 % um
sigma_vals = logspace(-4,-0.3, 192);
kD_vals = [0.3];                        % picoJ/um^2
kappa_vals = 1e-7;         % picoJ
epsilon_vals = -0.3e-3;              % picoJ/um^2
n0_vals = 1;                                    % fraction
alpha_i_vals = [1e-10, 1e-5, 1e-4, 1e-3];                            % fraction

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
                            if any(rr==rv) && any(ss==sv) ...
                                    && any(kk==kv) && any(ee==ev) ...
                                    && any(rr==rv) && any(pp==pv) ...
                                    && any(aa==av)
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



































