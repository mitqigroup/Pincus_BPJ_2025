% Get data

% location = 'task42/matlab.mat';
location = 'task25_QMQ/matlab.mat';
%location = 'isaac_old/task25/matlab.mat';
load(location)

ep_vals = -linspace(0.07,0.2,8)/2;
alpha_i0_vals = linspace(0.5,2,6)*1e-5;

[m,n] = ndgrid(ep_vals, alpha_i0_vals);
all_vals = [m(:),n(:)];

my_vals = all_vals(MyTaskID, :);

epsilon = my_vals(1)
alpha_i0 = my_vals(2)

% working so far:
% task17 
% task25
% task30
% task41
% task34 


newcolors = flip(["#003f5c", "#374c80", "#ef5675", "#bc5090", "#7a5195", ...
    "#ff764a", "#ffa600"]);

%% overall extraction
figure('Position',[488,173,782.5999999999999,589]);
hold on
rate_exp = [0.5, 1.5, 2.5, 2, 1.25]*1e3;
size_p = R_vals*2-18; % particle diameter accounting for 9 nm protein coating
rate = 1./(tw1+tw2);
% uptake = rate/max(rate(3:end))*max(rate_exp);
uptake = rate*2*3600*N_sites;
% plot(size_p(3:end), uptake(3:end), 'ro-', 'color', 'b', 'MarkerFaceColor','b',...
%     'Displayname', sprintf('$\\epsilon n_0 = %0.3g$', epsilon), ...
%     'MarkerSize', 12);
plot(size_p(5:end), uptake(5:end), 'ro-', 'color', 'b', 'MarkerFaceColor','b',...
    'Displayname', 'Our Model', ...
    'MarkerSize', 12);
% plot(size_p(5:end), uptake(5:end), 'ro-', 'color', 'b', 'MarkerFaceColor','b',...
%     'Displayname', 'Our Model', ...
%     'MarkerSize', 12);
size_p_exp = [14,30,50,74,100];
% rate_exp = [3, 4.5, 6, 4, 1.8];
% rate_exp = rate_exp/max(rate_exp)s;
plot(size_p_exp, rate_exp, 'rs:', 'color', 'r', 'MarkerFaceColor','r',...
    'DisplayName', 'Experiments','MarkerSize', 12);
xlabel('Particle diameter (nm)');
ylabel('Particles/cell at 2hrs');
% ylim([0,1.1])

% a1 = annotation('textbox', 'String',...
%     [sprintf('$R = %0.3g$ nm \n', R),...
%     sprintf('$2 \\sqrt{2 \\kappa/\\epsilon n_0} = %0.3g $ nm \n', ...
% 2*sqrt(2*kappa/-epsilon)),...
%     sprintf('$\\kappa = %0.3g$ $k_\\mathrm{B}T$ \n', kappa),...
%     sprintf('$d = %0.3g$ nm \n', d),...
%     sprintf('$\\tau = %0.3g$ s \n', 1/D),...
%     sprintf('$k_D = %0.3g$ $k_\\mathrm{B}T$/nm$^2$', kD)]);

a1 = annotation('textbox', 'String',...
    [sprintf('$\\epsilon n_0 = %0.3g $ mN/m \n',-epsilon/250*1000),...
    sprintf('$\\kappa = %0.3g$ $k_\\mathrm{B}T$ \n', kappa),...
    sprintf('$d = %0.3g$ $\\mu$m \n', d/1000),...
    sprintf('$\\tau = %0.3g$ s \n', tau),...
    sprintf('$k_D = %0.3g$ N/m \n', kD/250),...
    sprintf('$\\alpha_i(0) = %0.3g$', alpha_i0)]);
a1.Position = [0.161266666666667,0.456403269754769, ...
    0.379246153846154,0.218397820163489];

%%
opts = detectImportOptions("Lipowsky_uptake.csv");
data = readmatrix("Lipowsky_uptake.csv", opts);

p = 1-1e-1;
X = data(:,1);
Y = data(:,2)*(1/1.5e-3)*N_sites;
% Y = data(:,2);
pp = csaps(X,Y,p);
% pp = csaps(X,Y);
fnplt(pp);
% plot(data(:,1), data(:,2)*(1/1.5e-3))
% 
opts = detectImportOptions("Gao_data.csv");
data = readmatrix("Gao_data.csv", opts);

% p = 0.2;
X = data(:,1)*(sqrt(kappa/5e-3))*2-18;
Y = 2*3600./(data(:,2)*(kappa/(5e3*1e-2)))*N_sites;
% plot(X,Y, 'rx');

pp = csaps(X,Y,p);
% pp = csaps(X,Y);
fnplt(pp);
% 
axes1 = gca;
axes1.YScale = 'log';

legend(["Our Model", "Experiments [26]", "Model 1 [19]", "Model 2 [27]"])

xlim([10, 105])

%% different tau values
figure('Position',[488,173,782.5999999999999,589]);
hold on
rate_exp = [0.5, 1.5, 2.5, 2, 1.25]*1e3;
size_p = R_vals*2-18;
rate = 1./(tw1+tw2);
% uptake = rate/max(rate(3:end))*max(rate_exp);
tau_updated = [0.1,1,10,100];
for ii = 1:length(tau_updated)
    uptake = rate*2*3600*N_sites*tau/tau_updated(ii);
    plot(size_p(5:end), uptake(5:end), 'ro-', 'color', newcolors(ii),...
        'MarkerFaceColor',newcolors(ii),...
        'Displayname', sprintf('$\\tau = %0.3g$ s',tau_updated(ii)), ...
        'MarkerSize', 12);
end
size_p_exp = [14,30,50,74,100];
% rate_exp = [3, 4.5, 6, 4, 1.8];
% rate_exp = rate_exp/max(rate_exp)s;
plot(size_p_exp, rate_exp, 'rs:', 'color', 'r', 'MarkerFaceColor','r',...
    'DisplayName', 'Experiments','MarkerSize', 12);
xlabel('Particle Size (nm)');
ylabel('Particles/cell at 2hrs');
% ylim([0,1.1])

% 
% a1 = annotation('textbox', 'String',...
%     [sprintf('$\\epsilon n_0 = %0.3g $ mN/m \n',-epsilon/250*1000),...
%     sprintf('$\\kappa = %0.3g$ $k_\\mathrm{B}T$ \n', kappa),...
%     sprintf('$d = %0.3g$ $\\mu$m \n', d/1000),...
%     sprintf('$\\tau = %0.3g$ s \n', tau),...
%     sprintf('$k_D = %0.3g$ N/m \n', kD/250),...
%     sprintf('$\\alpha_i(0) = %0.3g$', alpha_i0)]);
a1.Position = [0.161266666666667,0.456403269754769, ...
    0.379246153846154,0.218397820163489];

axes1 = gca;
axes1.YScale = 'log';

legend
