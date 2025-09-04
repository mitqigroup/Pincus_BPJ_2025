clear variables

tic
if isfile('combined.mat')
    load('combined.mat')
else
    num_tasks = 96;
    E_all_toroid = [];
    alpha_A_vals_toroid = [];
    alpha_B_vals_toroid = [];
    S_A_vals_toroid = [];
    S_B_vals_toroid = [];
    Sigma_vals_toroid = [];
    rho_vals = [];
    
    E_all_linear = [];
    alpha_A_vals_linear = [];
    alpha_B_vals_linear = [];
    S_A_vals_linear = [];
    S_B_vals_linear = [];
    Sigma_vals_linear = [];
    
    E_all_nonlinear = [];
    alpha_A_vals_nonlinear = [];
    alpha_B_vals_nonlinear = [];
    S_A_vals_nonlinear = [];
    S_B_vals_nonlinear = [];
    Sigma_vals_nonlinear = [];

    for ii=0:num_tasks-1
        load_name = sprintf('data/task%i_results.mat', ii )
        data = load(load_name);
        load(load_name);
        E_all_toroid = cat(3,E_all_toroid, data.E_all_toroid_self);
        alpha_A_vals_toroid = horzcat(alpha_A_vals_toroid, data.alpha_A_vals_toroid_self);
        alpha_B_vals_toroid = horzcat(alpha_B_vals_toroid, data.alpha_B_vals_toroid_self);
        S_A_vals_toroid = horzcat(S_A_vals_toroid, data.S_A_vals_toroid_self);
        S_B_vals_toroid = horzcat(S_B_vals_toroid, data.S_B_vals_toroid_self);
        Sigma_vals_toroid = horzcat(Sigma_vals_toroid, data.Sigma_vals_toroid_self);
        rho_vals = horzcat(rho_vals, data.rho_vals_self);
        
        E_all_linear = cat(3,E_all_linear, data.E_all_linear_self);
        alpha_A_vals_linear = horzcat(alpha_A_vals_linear, data.alpha_A_vals_linear_self);
        alpha_B_vals_linear = horzcat(alpha_B_vals_linear, data.alpha_B_vals_linear_self);
        S_A_vals_linear = horzcat(S_A_vals_linear, data.S_A_vals_linear_self);
        S_B_vals_linear = horzcat(S_B_vals_linear, data.S_B_vals_linear_self);
        Sigma_vals_linear = horzcat(Sigma_vals_linear , data.Sigma_vals_linear_self);
        
        E_all_nonlinear = cat(3,E_all_nonlinear, data.E_all_nonlinear_self);
        alpha_A_vals_nonlinear = horzcat(alpha_A_vals_nonlinear, data.alpha_A_vals_nonlinear_self);
        alpha_B_vals_nonlinear = horzcat(alpha_B_vals_nonlinear, data.alpha_B_vals_nonlinear_self);
        S_A_vals_nonlinear = horzcat(S_A_vals_nonlinear, data.S_A_vals_nonlinear_self);
        S_B_vals_nonlinear = horzcat(S_B_vals_nonlinear, data.S_B_vals_nonlinear_self);
        Sigma_vals_nonlinear = horzcat(Sigma_vals_nonlinear, data.Sigma_vals_nonlinear_self);

    end
    
    parameter_set = data.parameter_set;
    save('combined.mat')
end
toc

% randomly shuffle the parameter set so that we have a more even distribution between
% each core on our node. We can unshuffle it later on if we want to. We MUST set the
% random number generator to the same value, since it must be identical between cores
rng('default');
rng(1);

size_params = size(parameter_set);

permutation_array = randperm(size_params(1));

% parameter_set_permuted = parameter_set;
% % parameter_set_unpermuted = parameter_set;
% 
% for ii=1:size_params(2)
%     parameter_set_permuted(:,ii) = parameter_set(permutation_array,ii);
% end
% parameter_set_original = parameter_set;
% parameter_set = parameter_set_permuted;

for jj = 1:length(permutation_array)
    kk = permutation_array(jj);
    parameter_set_unpermuted(kk,:) = parameter_set(jj,:);
    E_all_toroid_unpermuted(:,:,kk) = E_all_toroid(:,:,jj);
    alpha_A_vals_toroid_unpermuted(:,kk) = alpha_A_vals_toroid(:,jj);
    alpha_B_vals_toroid_unpermuted(:,kk) = alpha_B_vals_toroid(:,jj);
    S_A_vals_toroid_unpermuted(:,kk) = S_A_vals_toroid(:,jj);
    S_B_vals_toroid_unpermuted(:,kk) = S_B_vals_toroid(:,jj);
    Sigma_vals_toroid_unpermuted(:,kk) = Sigma_vals_toroid(:,jj);
    rho_vals_unpermuted(:,kk) = rho_vals(:,jj);
    
    E_all_linear_unpermuted(:,:,kk) = E_all_linear(:,:,jj);
    alpha_A_vals_linear_unpermuted(:,kk) = alpha_A_vals_linear(:,jj);
    alpha_B_vals_linear_unpermuted(:,kk) = alpha_B_vals_linear(:,jj);
    S_A_vals_linear_unpermuted(:,kk) = S_A_vals_linear(:,jj);
    S_B_vals_linear_unpermuted(:,kk) = S_B_vals_linear(:,jj);
    Sigma_vals_linear_unpermuted(:,kk) = Sigma_vals_linear(:,jj);
    
    E_all_nonlinear_unpermuted(:,:,kk) = E_all_nonlinear(:,:,jj);
    alpha_A_vals_nonlinear_unpermuted(:,kk) = alpha_A_vals_nonlinear(:,jj);
    alpha_B_vals_nonlinear_unpermuted(:,kk) = alpha_B_vals_nonlinear(:,jj);
    S_A_vals_nonlinear_unpermuted(:,kk) = S_A_vals_nonlinear(:,jj);
    S_B_vals_nonlinear_unpermuted(:,kk) = S_B_vals_nonlinear(:,jj);
    Sigma_vals_nonlinear_unpermuted(:,kk) = Sigma_vals_nonlinear(:,jj);
end

parameter_set = parameter_set_unpermuted;
E_all_toroid = E_all_toroid_unpermuted;
alpha_A_vals_toroid = alpha_A_vals_toroid_unpermuted;
alpha_B_vals_toroid = alpha_B_vals_toroid_unpermuted;
S_A_vals_toroid = S_A_vals_toroid_unpermuted;
S_B_vals_toroid = S_B_vals_toroid_unpermuted;
Sigma_vals_toroid = Sigma_vals_toroid_unpermuted;
rho_vals = rho_vals_unpermuted;

E_all_linear = E_all_linear_unpermuted;
alpha_A_vals_linear = alpha_A_vals_linear_unpermuted;
alpha_B_vals_linear = alpha_B_vals_linear_unpermuted;
S_A_vals_linear = S_A_vals_linear_unpermuted;
S_B_vals_linear = S_B_vals_linear_unpermuted;
Sigma_vals_linear = Sigma_vals_linear_unpermuted;

E_all_nonlinear = E_all_nonlinear_unpermuted;
alpha_A_vals_nonlinear = alpha_A_vals_nonlinear_unpermuted;
alpha_B_vals_nonlinear = alpha_B_vals_nonlinear_unpermuted;
S_A_vals_nonlinear = S_A_vals_nonlinear_unpermuted;
S_B_vals_nonlinear = S_B_vals_nonlinear_unpermuted;
Sigma_vals_nonlinear = Sigma_vals_nonlinear_unpermuted;

%% base case
% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = [0.05,0.5];                 % um
sigma_vals = logspace(-6,-2,12);                % surface fraction
%d_vals = sqrt(R_vals.^2./sigma_vals);           % um
kD_vals = 300/10^12*1e9;                        % picoJ/um^2
kD_base = 300/10^12*1e9;                        % picoJ/um^2
zeta_vals = logspace(-6,-3,12);                  % dimensionless
epsilon_vals = -zeta_vals*kD_base;              % picoJ/um^2
n0_vals = 1;                                    % fraction
kappa_vals = 1e-19*1e12;         % picoJ
% kappa_vals = logspace(-21,-15, 60)*1e12;
alpha_i_vals = [1e-6,1e-5,1e-4,1e-3,1e-2];                            % fraction
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
%                             if ee==12 && aa==1 && rr==1 && ss == 1
%                                 slice = [slice, ii];
%                             end
%                             if ee==6 && aa==1 && rr==1 && ss == 1
%                                 slice = [slice, ii];
%                             end
%                             if ee==12 && aa==2 && rr==2 && ss == 7
%                                 slice = [slice, ii];
%                             end
                            if ee==10 && aa==2 && rr==1 && ss==2
                                slice = [slice, ii];
                            end
%                             if ee==12 && aa==2 && rr==1 && ss==6
%                                 slice = [slice, ii];
%                             end
%                             if aa==1 && rr==2 && ss == 1
%                                 slice = [slice, ii];
%                             end
                        end
                    end
                end
            end
        end
    end
end
size_params = size(parameter_set);
base_index = slice(1);
newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];


% epsilon = parameter_set(base_index,1);
% n0 = parameter_set(base_index,2);
% d = parameter_set(base_index,3);
% R = parameter_set(base_index,4);
% kD = parameter_set(base_index,5);
% kappa = parameter_set(base_index,6);
% alpha_i = parameter_set(base_index,7);
% zeta = epsilon*n0/kD;
% sigma = pi*R^2/d^2;
% R_vals = parameter_set(slice,4)';
% d_vals = parameter_set(slice,3)';
% kappa_vals = parameter_set(slice, 6)';
% zeta_vals = parameter_set(slice,1).*parameter_set(slice,2)./parameter_set(slice,5);
% 
% changing_value = 3;
% param_1_vals = pi*R^2./parameter_set(slice,changing_value).^2;
% % param_1_vals = parameter_set(slice,changing_value);
% E_all = E_all_total(:, slice);
% alpha_A_vals = alpha_A_vals_total(slice);
% alpha_B_vals = alpha_B_vals_total(slice);
% phi_vals = phi_vals_total(slice);
% % h_phi_vals = h_phi_vals_total(slice);
% S_A_vals = S_A_vals_total(slice);
% S_B_vals = S_B_vals_total(slice);
% Sigma_vals = Sigma_vals_total(slice);
% sigma_vals = pi*R_vals.^2./d_vals.^2;
% lambda_vals = sqrt(kappa_vals./(kD*alpha_B_vals));
% 
% solution_set = solution_set_total(slice);
% 
% values_structure.E_all = E_all_total(:, slice);
% values_structure.alpha_A_vals = alpha_A_vals_total(slice);
% values_structure.alpha_B_vals = alpha_B_vals_total(slice);
% values_structure.phi_vals = phi_vals_total(slice);
% % values_structure.h_phi_vals = h_phi_vals_total(slice);
% values_structure.S_A_vals = S_A_vals_total(slice);
% values_structure.S_B_vals = S_B_vals_total(slice);
% values_structure.Sigma_vals = Sigma_vals_total(slice);
% values_structure.sigma_vals = pi*R_vals.^2./d_vals.^2;
% values_structure.lambda_vals = sqrt(kappa_vals./(kD*alpha_B_vals));
% 
% values_structure.R_vals = parameter_set(slice,4)';
% values_structure.d_vals = parameter_set(slice,3)';
% values_structure.kappa_vals = parameter_set(slice, 6)';
% values_structure.zeta_vals = parameter_set(slice,1).*parameter_set(slice,2)./parameter_set(slice,5);
% values_structure.param_1_vals = param_1_vals;
% 
% values_structure.epsilon = parameter_set(base_index,1);
% values_structure.n0 = parameter_set(base_index,2);
% values_structure.d = parameter_set(base_index,3);
% values_structure.R = parameter_set(base_index,4);
% values_structure.kD = parameter_set(base_index,5);
% values_structure.kappa = parameter_set(base_index,6);
% values_structure.alpha_i = parameter_set(base_index,7);
% values_structure.zeta = epsilon*n0/kD;
% values_structure.sigma = R^2/d^2;

colours_cb = ["#27E0D8", "#648FFF", "#1D0C6B", "#DC267F", "#FE6100", ...
    "#FFB000", "#016D24"];

% close all

% xaxis_name = '$\alpha_i$';
xaxis_name = '$\sigma$';
% xaxis_name = '$R$ ($\mu$m)';
% xaxis_name = '$\kappa$ (pJ)';

% xscale = 'linear';
xscale = 'log';



%%%
f1 = figure('Position',[400,100,800,600]);
hold on
% xlim([0,180])
xlabel('$\phi$')
ylabel('$\Delta E$')
% ylabel('$\Delta E/\mathrm{max}(|\Delta E|)$')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
% jj = 16;
for jj = slice
epsilon = parameter_set(jj,1)
n0 = parameter_set(jj,2);
d = parameter_set(jj,3);
R = parameter_set(jj,4);
kD = parameter_set(jj,5);
kappa = parameter_set(jj,6);
alpha_i = parameter_set(jj,7);
sigma = pi*R^2/d^2;
zeta = -epsilon*n0/kD;
anno_string = [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
    sprintf('$\\kappa = %0.2g$ pJ \n', kappa),...
    sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
    sprintf('$R = %0.2g$ $\\mu$m \n', R),...
    sprintf('$\\zeta = %0.2g$ \n', zeta),...
    sprintf('$\\sigma = %0.2g$ \n', sigma)];
E_unwrapped = [kD/2*alpha_i^2*d^2/(1+alpha_i),0,0,kD/2*alpha_i^2*d^2/(1+alpha_i),0,0,nan];
% E_unwrapped_toroid = zeros(size(E_all_toroid));
% E_unwrapped_toroid([1,4],:) = repmat(kD/2*alpha_B_vals_toroid.^2.*...
%     (d^2-R^2*sin(phi_vals).^2)./(1+alpha_B_vals_toroid.^2),2,1);
% E_unwrapped_linear = zeros(size(E_all_linear));
% E_unwrapped_linear([1,4],:) = repmat(kD/2*alpha_B_vals_linear.^2.*...
%     (d^2-R^2*sin(phi_vals).^2)./(1+alpha_B_vals_linear.^2),2,1);
% ab_end = alpha_B_vals_toroid(1);
z_vals = (1-cos(phi_vals));
% for ii=[1]
for ii=1:6
%     plot(rad2deg(phi_vals), E_all(ii,:)-E_all(ii,end), ...
%         strcat(colours(ii),lines(ii)))
%     plot(rad2deg(phi_vals), E_all_linear(ii,:), ...
%         strcat(colours(ii),'-'))
%     plot(rad2deg(phi_vals), E_all_linear(ii,:,jj)-E_unwrapped(ii),...
%         '-','displayname', 'linear')
%     plot(rad2deg(phi_vals), E_all_toroid(ii,:,jj)-E_unwrapped(ii),...
%         '--','displayname', 'toroid')
%     plot(rad2deg(phi_vals), E_all_nonlinear(ii,:,jj)-E_unwrapped(ii),...
%         '-','displayname', 'nonlinear')
    E_n = E_all_nonlinear(ii,:,jj);
%     if ii==1||ii==4||ii==6
%         E_t = E_all_toroid(ii,:,jj);
%         E_n(E_n>E_t) = nan;
%     end
    E_n = E_n - E_unwrapped(ii);
    p1 = plot(rad2deg(phi_vals), E_n,...
        '-','displayname', 'nonlinear');
%     p1 = plot(rad2deg(phi_vals), E_n/max(abs(E_n)),...
%         '-','displayname', sprintf('$\\sigma = %0.3g$', sigma));
%     p1 = plot(rad2deg(phi_vals), E_n/max(abs(E_n)),...
%         '-','displayname', sprintf('$\\zeta = %0.3g$', zeta));
    if ii==1
        p1.LineWidth = 3;
    end
    if ii==6
%         p1 = plot(rad2deg(phi_vals), E_all_nonlinear(2,:,jj)+ E_all_nonlinear(5,:,jj),...
%         '-','displayname', 'nonlinear');
        p1 = plot(rad2deg(phi_vals), E_all_nonlinear(2,:,jj)+ E_all_nonlinear(4,:,jj)+...
            E_all_nonlinear(5,:,jj) - E_unwrapped(1),...
        '-','displayname', 'nonlinear');
        x = phi_vals(1:3:end);
        ab = (1+sigma*(2*(1-cos(x))-sin(x).^2))*(1+alpha_i)-1;
        ytest = kD/2*(d^2-pi*R^2*sin(x).^2).*ab.^2./(1+ab)-kD/2*d^2*alpha_i^2/(1+alpha_i);
%         constant = mean(S_A_vals_nonlinear(:,jj))./(1-cos(x'));
        p1 = plot(rad2deg(x), 2*pi*R^2*epsilon*n0*(1-cos(x')), 'k^', 'MarkerSize',6);
        p1 = plot(rad2deg(x), 4*pi*kappa*(1-cos(x')), 'ko', 'MarkerSize',6);
        p1 = plot(rad2deg(x), ytest, 'ks', 'MarkerSize',6);
    end
%     plot(rad2deg(phi_vals), E_n/max(abs(E_n)),...
%         '-','displayname', 'nonlinear')
%     plot(rad2deg(phi_vals), E_all_nonlinear(ii,:,jj)-E_unwrapped(ii),...
%         '-','displayname', sprintf('$\\sigma = %0.3g$', sigma))
%     plot(z_vals(1:end-1), diff((E_all_linear(ii,:)-E_unwrapped(ii)))./diff(z_vals),...
%         '-','displayname', sprintf('$\\epsilon n_0 R^2/(2 \\kappa) = %0.2g$, $\\alpha_i = %0.2g$',...
%         ep_bar, alpha_i))
%     plot(rad2deg(phi_vals), E_all_linear(ii,:)-E_unwrapped(ii),...
%         '-','displayname', sprintf('$\\zeta/ \\alpha_B = %0.2g$, linear', zeta/ab_end))
%     plot(rad2deg(phi_vals), E_all_linear(ii,:)-E_unwrapped_linear(ii,:),...
%         '-','displayname', sprintf('$\\zeta/ \\alpha_B = %0.2g$, linear', zeta/ab_end))
%     plot(rad2deg(phi_vals), E_all_toroid(ii,:)-E_unwrapped_toroid(ii,:), ...
%         '--','displayname', sprintf('$\\zeta/ \\alpha_B = %0.2g$, toroid', zeta/ab_end))
%     plot(rad2deg(phi_vals), E_all_toroid(ii,:)-E_unwrapped(ii), ...
%         '-','displayname', sprintf('$\\sigma = %g$', sigma))
%     plot(rad2deg(phi_vals), E_all_toroid(ii,:)-E_unwrapped(ii), ...
%         '--','displayname', sprintf('$\\epsilon/\\Sigma = %g$, toroid', -epsilon/(kD*alpha_i)))

end
end
colororder(newcolors)

xlim([0,180])
annotation('textbox', 'String', anno_string)

% legend
legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
    '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$',...
    '$E_\mathrm{adhesion}+E_\mathrm{bend,A}$',...
    '$2 \pi R^2 \epsilon n_0 (1-\cos \phi)$', ...
    '$4 \pi \kappa (1-\cos \phi)$'}, 'Box','off',...
    'location', 'best')

%% analyse data from each, we want local minima and maxima and end points.
% for jj = 1:length(parameter_set)
counter = 0;
for jj = slice
    counter = counter+1;
    d = parameter_set(jj,3);
    kD = parameter_set(jj,5);
    alpha_i = parameter_set(jj,7);
    E = squeeze(E_all_nonlinear(1,:,jj))-kD/2*alpha_i^2*d^2/(1+alpha_i);
    diff_E = diff(E);
    % find points where diff changes sign, these are stationary points.
    locs_pos = [];
    locs_neg = [];
    for ii = 2:length(diff_E)
        if diff_E(ii)<0&&diff_E(ii-1)>0
            locs_pos = [locs_pos, ii];
        elseif diff_E(ii)>0&&diff_E(ii-1)<0
            locs_neg = [locs_neg, ii];
        end
    end
    [min_val(counter), min_loc(counter)] = min(E);
    [max_val, max_loc] = max(E);
    E_low_phi = E(end);
    E_high_phi = E(1);

    % find a few points near the minimum and do a quadratic fit, if the
    % minimum isn't at the boundaries
    if min_loc(counter)~=1&&min_loc(counter)~=length(phi_vals)
        % get some points either side
        min_range =(min_loc(counter)-1):(min_loc(counter)+1);
        phi_near_min = phi_vals(min_range);
        E_near_min = E(min_range);
        pol_coeff = polyfit(phi_near_min, E_near_min, 2);
%         % graph if we want
%         figure();
%         hold on
% %         plot(phi_near_min, E_near_min)
%         plot(phi_vals, E)
%         x = linspace(phi_near_min(1),phi_near_min(end));
%         plot(x, polyval(pol_coeff, x))
        % make sure that the second derivative isnt negative
        if pol_coeff(1)<0
            error('second derivative is negative!');
        end
        % then get the minima
        phi_min(counter) = -pol_coeff(2)/(2*pol_coeff(1));


    else
        phi_min(counter) = phi_vals(min_loc(counter));
    end
    

    out = free_shape_nonlinear_free_h(R,d,phi_min(counter),...
        kappa, alpha_B_vals_nonlinear(min_loc(counter), jj),1);
    r = out.y(2,:);
    h = out.y(3,:);
    figure('Position',[400,100,800,600]);
    hold on
    axis equal
    xlabel('$r$')
    ylabel('$h$')
    p1 = plot(r, h-h(1), 'r-', 'displayname', 'free surface');
    p1.Color = newcolors(color_number);
    t = linspace(0,2*pi,1000);
    x = cos(t)*R;
    % y = sin(t)*R+(R*cos(phi)+h(1));
    y = sin(t)*R+R*cos(phi_min(counter));
    p1 = plot(x,y,'r:', 'displayname', 'microbead');
    p1.Color = newcolors(color_number);

    t = linspace(-pi/2,-pi/2+phi_min(counter),1000);
    x = cos(t)*R;
    % y = sin(t)*R+(R*cos(phi)+h(1));
    y = sin(t)*R+R*cos(phi_min(counter));
    p1 = plot(x,y,'r-', 'displayname', 'attached surface');
    p1.Color = newcolors(color_number);

    axis off

end
d_vals = parameter_set(slice,3);
R_vals = parameter_set(slice,4);
epsilon = parameter_set(slice(1),1);
n0 = parameter_set(slice(1),2);
d = parameter_set(slice(1),3);
R = parameter_set(slice(1),4);
kD = parameter_set(slice(1),5);
kappa = parameter_set(slice(1),6);
alpha_i = parameter_set(slice(1),7);
sigma_vals = pi*R_vals.^2./d_vals.^2;
% anno_string = [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
%     sprintf('$\\kappa = %0.2g$ pJ \n', kappa),...
%     sprintf('$\\alpha_i = %0.2g$ \n', alpha_i)];
anno_string = [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
    sprintf('$\\kappa = %0.2g$ pJ \n', kappa),...
    sprintf('$\\epsilon n_0 = %0.2g$ \n', epsilon*n0),...
    sprintf('$R = %0.2g$ pJ \n', R)];
% figure()
% axes1 = gca;
% hold on
% xlabel('$\sigma$')
% ylabel('$\phi$')
% axes1.XScale = 'log';
% yticks([0, pi/4, pi/2, 3*pi/4, pi])
% yticklabels({'$0$', '$\pi/4$', '$\pi/2$', '$3 \pi/4$', '$\pi$'})
% ylim([0,pi+0.01])
% % plot(sigma_vals, phi_vals(min_loc),...
% %     'displayname', sprintf('$\\epsilon n_0 = %0.3g$', epsilon*n0))
% % plot(sigma_vals, phi_min,...
% %     'displayname', sprintf('$\\epsilon n_0 = %0.3g$', epsilon*n0))
% % plot(sigma_vals, phi_min, '-',...
% %     'displayname', sprintf('$\\epsilon n_0 = %0.3g$, $R = %0.3g$ $\\mu$m', epsilon*n0, R))
% plot(sigma_vals, phi_min, '-',...
%     'displayname', sprintf('$\\alpha_i = %0.3g$', alpha_i))
% annotation('textbox', 'String', anno_string)
% legend

%%
% newcolors = ["#FFB000","#27E0D8","#1D0C6B","#DC267F","#FE6100","#648FFF","#016D24"];
% newcolors = ["#016D24","#648FFF","#FE6100","#DC267F","#1D0C6B","#27E0D8","#FFB000"];
newcolors = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];
% newcolors = flip(["#016D24","#FE6100","#DC267F"]);
% figure();
% plot(1:7,ones(1,7)'*[1:7])
colororder(newcolors)

