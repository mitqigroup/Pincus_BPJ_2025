clear variables
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');

% lookup_table = load("lookup_table.mat")
% lookup_table = load("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Nonlinear_solutions/gridded_interp.mat");
lookup_table = load("gridded_interp.mat");

% all lengths in um, all energies in pJ

% kB = 1.380649e-23;
kB = 1.3e-23;
NA = 6.0221408e+23;
T = 298.15; % K
% R = 0.1; % um
kD = 0.3; %pJ/um^2
% kD = 10e-5; %pJ/um^2
% epsilon = -1e-7;
kappa = 25*kB*T*1e12; % pJ
% mu = 8e-6; % pJ/um^2
mu = 0; % pJ/um^2
initial_tension = 1e-7; % pJ/um^2
A_total = 140; % um^2
phi_vals = deg2rad(linspace(1,179,200));

Sigma = initial_tension;

% epsilon_vals = -logspace(-8,-4,10);
% R_vals = [0.025,0.05,0.1,0.25];
% N_vals = [1,2,5,10,30,50,100,200];
N_vals = [1];
epsilon_vals = -3e-5; % pJ/um^2
R_vals = 0.1;
% N_vals = 2;
% E = zeros(length(phi_vals));
% Sig = zeros(size(phi_vals));
% Ap = zeros(size(phi_vals));

sigma_vals = N_vals.*R_vals.^2*pi/A_total

for ee=1:length(epsilon_vals)
epsilon = epsilon_vals(ee);
for rr=1:length(R_vals)
R = R_vals(rr);
tic
for nn=1:length(N_vals)
N_particles = N_vals(nn);
Ai = A_total/N_particles;
A = Ai;
d = sqrt(A);
alpha_i = get_alpha(initial_tension, A_total*1e-12, kappa*1e-12, kD, kB*T);
Sigma = initial_tension;
for pp=1:length(phi_vals)
    phi = phi_vals(pp);
    % find the correct value of Sigma by iteration
    counter = 0;
    while true
        counter = counter+1;
        Sigma_prev = Sigma;
        sig_bar = Sigma*R^2/kappa;
        % use interpolant
        delA = lookup_table.F_delA(sig_bar,phi,d/R)*R^2;
        % full solution
%         out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, false);
        A = Ai+pi*R^2*(1-cos(phi))+delA;
        alpha = (A/Ai*(1+alpha_i)-1);
        a = kB*T/(8*pi*kappa*1e-12);
        b = 0.1*A/kappa;
        c = kD;
        Sigma = (a*b*c*lambertw(exp(alpha/a+1/(a*b*c))/(a*b*c))-1)/b;
%         Sigma = fzero(@(x) get_alpha(x, A*1e-12, kappa*1e-12, kD, kB*T)-alpha, Sigma*10);
        if abs(Sigma_prev-Sigma)/Sigma<1e-6
            break
        elseif counter>20
            fprintf('over 20 iterations, rel err is %d \n', abs(Sigma_prev-Sigma)/Sigma);
            break
        end
    end
    out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, false);
    hphi(ee,rr,nn,pp) = out.y(3,end)-out.y(3,1)+R*(1-cos(phi));
    % calculate values
    sig_bar = Sigma*R^2/kappa;
    delA = lookup_table.F_delA(Sigma*R^2/kappa,phi,d/R)*R^2;
    A = Ai+pi*R^2*(1-cos(phi))+delA;
    A_cap = 2*pi*R^2*(1-cos(phi));
    alpha = (A/Ai*(1+alpha_i)-1);
    E_bend_free = lookup_table.F_Ebend(Sigma*R^2/kappa,phi,d/R)*pi*kappa;
    E_shear = 2*pi*mu*R^2*2*(-1/4*sin(phi/2).^4-1/2*sin(phi/2).^2-log(cos(phi/2))) ...
        + mu*(A_cap-pi*R^2*sin(phi).^2).*log(A_cap./(pi*R^2*sin(phi).^2));
    E(ee,rr,nn,pp) = epsilon*(1-cos(phi))+Sigma*(A-Ai)+E_bend_free+4*pi*kappa*(1-cos(phi))+E_shear;
    Sig(ee,rr,nn,pp) = Sigma;
    Ap(ee,rr,nn,pp) = A;
end
end
toc
end
end

%%
figure();
plot(phi_vals, squeeze(E));

figure();
plot(phi_vals, squeeze(Sig));

% %%
% figure('Position',[201,158,770,579]);
% hold on
% for N_index = 1:length(N_vals)
% E_plot = squeeze(E(:,1,N_index,:));
% Sig_plot = squeeze(Sig(:,1,N_index,:));
% hphi_plot = squeeze(hphi(:,1,N_index,:));
% movmedian_range = 11;
% movmedian_threshold = 0.1;
% N_spline_points = 10000;
% use_weights = false;
% outliers = true(size(E_plot));
% for ii=1:2
%     [outliers,L,U,C] = isoutlier(E_plot, 'movmedian',movmedian_range, ...
%         'ThresholdFactor',movmedian_threshold);
%     outliers(1:floor(movmedian_range/2)) = false;
%     outliers(end-floor(movmedian_range/2):end) = false;
%     E_plot(outliers) = nan;
% end
% xx = linspace(0,pi,N_spline_points);
% E_spline = fit_spline(phi_vals, E_plot, xx, false(size(E_plot)));
% 
% % get minima in phi
% % [local_minima, min_promenance] = islocalmin(E_spline);
% % phi_local_min(N_index) = min(xx(local_minima));
% % yy = fit_spline(phi_vals, Sig_plot, xx, false(size(Sig_plot)));
% % Sig_local_min(N_index) = min(yy(local_minima));
% % E_local_min(N_index) = min(E_spline(local_minima));
% % yy = fit_spline(phi_vals, hphi_plot, xx, false(size(hphi_plot)));
% % hphi_local_min(N_index) = min(yy(local_minima));
% 
% plot(xx/pi, E_spline, 'LineWidth',2)
% % plot(phi_vals/pi, E_plot, 'LineWidth',2)
% plot(phi_vals/pi, squeeze(E(:,1,N_index,:)), '--','LineWidth',0.5)
% % plot(phi_vals(~outliers)/pi, E_plot(~outliers))
% end
% 
% figure();
% hold on
% xlabel('$N_\mathrm{particles}$')
% ylabel('$\phi/\pi$')
% plot(N_vals, phi_local_min/pi);
% yyaxis right
% ylabel('$\alpha_B$ ($\mu$N/m)')
% plot(N_vals, Sig_local_min*1e6);

%% functions

% get_alpha(1e-3, 8600*1e-12, 1e-19, 0.3, 1.380649e-23)

function alpha = get_alpha(Sigma, A, kappa, kD, kBT)
    % all SI units
    % Evans 2000 Elasticity Lipid Bilayers
    alpha = kBT/(8*pi*kappa)*log(1+0.1*A*Sigma/kappa)+Sigma/kD;
    % just kD
%     alpha = Sigma/kD;
end