% combine data from several runs

MyTaskID = 0;
NumberOfTasks = 1;

% check that the environment variables have been read in correctly
if ~(exist('MyTaskID', 'var')&&exist('NumberOfTasks', 'var'))
    error('Environment variables not set correctly')
end

load_name = sprintf('data/task%i_results.mat', MyTaskID);
load(load_name);

E_all_toroid = nan(6,length(phi_vals));
alpha_A_vals_toroid = nan(1,length(phi_vals));
alpha_B_vals_toroid = nan(size(alpha_A_vals_toroid));
S_A_vals_toroid = nan(size(alpha_A_vals_toroid));
S_B_vals_toroid = nan(size(alpha_A_vals_toroid));
Sigma_vals_toroid = nan(size(alpha_A_vals_toroid));
rho_vals = nan(size(alpha_A_vals_toroid));

E_all_linear = nan(6,length(phi_vals));
alpha_A_vals_linear = nan(1,length(phi_vals));
alpha_B_vals_linear = nan(size(alpha_A_vals_toroid));
S_A_vals_linear = nan(size(alpha_A_vals_toroid));
S_B_vals_linear = nan(size(alpha_A_vals_toroid));
Sigma_vals_linear = nan(size(alpha_A_vals_toroid));

E_all_nonlinear = nan(6,length(phi_vals));
alpha_A_vals_nonlinear = nan(1,length(phi_vals));
alpha_B_vals_nonlinear = nan(size(alpha_A_vals_toroid));
S_A_vals_nonlinear = nan(size(alpha_A_vals_toroid));
S_B_vals_nonlinear = nan(size(alpha_A_vals_toroid));
Sigma_vals_nonlinear = nan(size(alpha_A_vals_toroid));

% This is serial so only run if MyTaskID = 0
if MyTaskID==0
    for taskID = 0:NumberOfTasks-1
        load_name = sprintf('data/task%i_results.mat', taskID);
        load(load_name);

        my_slice = my_set_min:my_set_max;

        E_all_toroid(:,my_slice) = E_all_toroid_self;
        alpha_A_vals_toroid(my_slice) = alpha_A_vals_toroid_self;
        alpha_B_vals_toroid(my_slice) = alpha_B_vals_toroid_self;
        S_A_vals_toroid(my_slice) = S_A_vals_toroid_self;
        S_B_vals_toroid(my_slice) = S_B_vals_toroid_self;
        Sigma_vals_toroid(my_slice) = Sigma_vals_toroid_self;
        rho_vals(my_slice) = rho_vals_self;
        
        E_all_linear(:,my_slice) = E_all_linear_self;
        alpha_A_vals_linear(my_slice) = alpha_A_vals_linear_self;
        alpha_B_vals_linear(my_slice) = alpha_B_vals_linear_self;
        S_A_vals_linear(my_slice) = S_A_vals_linear_self;
        S_B_vals_linear(my_slice) = S_B_vals_linear_self;
        Sigma_vals_linear(my_slice) = Sigma_vals_linear_self;
        
        E_all_nonlinear(:,my_slice) = E_all_nonlinear_self;
        alpha_A_vals_nonlinear(my_slice) = alpha_A_vals_nonlinear_self;
        alpha_B_vals_nonlinear(my_slice) = alpha_B_vals_nonlinear_self;
        S_A_vals_nonlinear(my_slice) = S_A_vals_nonlinear_self;
        S_B_vals_nonlinear(my_slice) = S_B_vals_nonlinear_self;
        Sigma_vals_nonlinear(my_slice) = Sigma_vals_nonlinear_self;

    end
end

f1 = figure('Position',[400,100,700,500]);
hold on
% xlim([0,180])
xlabel('$\phi$')
ylabel('$\Delta E$')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
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
%     plot(rad2deg(phi_vals), E_all_linear(ii,:)-E_unwrapped(ii),...
%         '-','displayname', 'linear')
%     plot(rad2deg(phi_vals), E_all_toroid(ii,:)-E_unwrapped(ii),...
%         '--','displayname', 'toroid')
    plot(rad2deg(phi_vals), E_all_nonlinear(ii,:)-E_unwrapped(ii),...
        '-.','displayname', 'nonlinear')
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

legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
    '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
    'location', 'best')

%%
load("data/task0_results.mat")

figure('Position',[400,100,800,600]);
hold on
axis equal
xlabel('$r$')
ylabel('$h$')
plot(r_nonlin, h_nonlin-h_nonlin(1), 'b-', 'displayname', 'free surface');
t = linspace(-pi/2,pi/2,1000);
x = cos(t)*R;
% y = sin(t)*R+(R*cos(phi)+h(1));
y = sin(t)*R+R*cos(phi);
plot(x,y,'b:', 'displayname', 'microbead')
xlim([0,5*R])
