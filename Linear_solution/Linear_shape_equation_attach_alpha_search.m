clear variables

% constants
R = 0.3;                  % um
sigma = 0.01;          % surface fraction
d = sqrt(R^2/sigma);    % um
% phi = pi/12;
kD  = 300/10^12*1e9;    % picoJ/um^2
zeta = 0.02;            % dimensionless
epsilon = -zeta*kD;     % picoJ/um^2 
% epsilon = -1;
n0 = 1;                 % fraction
% kappa = 1e-19*1e12;     % picoJ
kappa = 1e-17*1e12;     % picoJ
alpha_i = -0.01;        % fraction
zeta = epsilon*n0/kD;   % dimensionless
N = 3e3;                % number of points in quadrature

%% just plot the constraint

alpha_A = 0.01;
alpha_B = 0.01;
phi = deg2rad(25);
h_phi = -0.11;

% plot constraints
const = [epsilon, n0, d, R, kD, kappa, alpha_i, N];
inp = [alpha_A, alpha_B, phi, h_phi];
[c,ceq] = lipid_con_bend(inp, const);

alpha_A_vals = linspace(-0.1, 0.1, 20);
alpha_B_vals = linspace(-0.01, max(alpha_i*2,0.01), 50);

for ii=1:length(alpha_A_vals)
    alpha_A = alpha_A_vals(ii);
    for jj=1:length(alpha_B_vals)
        alpha_B = alpha_B_vals(jj);
        inp = [alpha_A, alpha_B, phi, h_phi];
        [~,ceq(ii,jj)] = lipid_con_bend(inp, const);
        E(ii,jj) = stretch_bend_min(inp, const);
    end
end

%%

figure();
hold on
xlabel('$\alpha_A$')
ylabel('$\alpha_B$')
% surf(alpha_A_vals, alpha_B_vals, ceq')
contour(alpha_A_vals, alpha_B_vals, ceq',[0,0],...
    '-r', 'linewidth', 1.5, 'DisplayName', 'Constraint')
contour(alpha_A_vals, alpha_B_vals, E', linspace(-1e-3, -1e-5, 5),...
    '-b', 'linewidth', 1.5, 'DisplayName', 'Objective')

legend('Location','northeast')

figure();
hold on
annotation('textbox', [0.3,0.6,0.4,0.2], 'String',...
    [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
    sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
    sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
    sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
    sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta),...
    sprintf('$\\kappa = %0.2g k_\\mathrm{D} R^2$ \n', kappa/(kD*R^2)),...
    sprintf('$\\phi = %0.2g^{\\circ}$ \n', rad2deg(phi)),...
    sprintf('$h_\\phi = %0.2g$ $\\mu$m \n', h_phi)],...
    'interpreter', 'latex','FontSize',16,...
    'FitBoxToText','on','LineStyle','none', 'Color','b')

figure();
hold on
xlabel('$\alpha_A$')
ylabel('$\alpha_B$')
surf(alpha_A_vals, alpha_B_vals, E', 'linewidth', 1.5)

% figure();
% hold on
% xlabel('$\alpha_B$')
% plot(alpha_B_vals, ceq(1,:))
% plot([min(alpha_B_vals), max(alpha_B_vals)], [0,0], 'k--')

% for ii = 1:length(alpha_A_vals)
%     alpha_A = alpha_A_vals(ii);
%     inp = [alpha_A, 0.1, phi, h_phi];
%     alpha_B_vals(ii) = fzero(@(x) get_constraint(x, inp, const), -0.1);
% end
% 
% figure();
% hold on
% plot(alpha_A_vals, alpha_B_vals)

function f = stretch_bend_min(y, const)
    % function of free energy to minimise

    epsilon = const(1);
    n0 = const(2);
    d = const(3);
    R = const(4);
    kD = const(5);
    kappa = const(6);
    alpha_i = const(7);
    N = const(8);
    
    alpha_A = y(1);
    alpha_B = y(2);
    phi = y(3);
    h_phi = y(4);
    
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
    r = linspace(r_phi, d/2,N);
    
    [~,~,S_B, ~, lap_h, ~] = free_shape_linear_fixed_h(...
        r, r_phi, d, phi, kappa, Sigma, h_phi);
    
%     hderiv = C(1)./r+C(3)/lambda*besseli(1,r/lambda)-C(4)/lambda*besselk(1,r/lambda);
%     lap_h = C(3)/lambda^2*besseli(0,r/lambda)+C(4)/lambda^2*besselk(0,r/lambda);
    S_A = 2*pi*R^2*(1-cos(phi));
    
    % stretching, adhesion and bending energy
    f = epsilon*n0*S_A./(1+alpha_A) ...
      + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
      + kappa/2*2*pi*trapz(r, r.*lap_h.^2) + 4*pi*kappa*(1-cos(phi));
% 
%     f = zeta*S_A./(1+alpha_A) ...
%       + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B));

%     f = epsilon*n0*S_A./(1+alpha_A) ...
%       + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B));

end

function [c,ceq] = lipid_con_bend(y, const)

    epsilon = const(1);
    n0 = const(2);
    d = const(3);
    R = const(4);
    kD = const(5);
    kappa = const(6);
    alpha_i = const(7);
    N = const(8);
    
    alpha_A = y(1);
    alpha_B = y(2);
    phi = y(3);
    h_phi = y(4);
    
    Sigma = kD*alpha_B;
    r_phi = sin(phi)*R;
    
    r = linspace(r_phi, d/2,N);
    [~,~,S_B, ~, ~, ~] = free_shape_linear_fixed_h(...
        r, r_phi, d, phi, kappa, Sigma, h_phi);
    
    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(phi));

    c(1) = 0;

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end

function f = get_constraint(alpha_B, y, const)
    
    [~, f] = lipid_con_bend([y(1), alpha_B, y(3), y(4)],...
        const);

end