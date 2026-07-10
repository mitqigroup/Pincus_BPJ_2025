clear variables

% constants
R = 0.05;                  % um
sigma = 0.01;          % surface fraction
d = sqrt(R^2/sigma);    % um
% phi = pi/12;
kD  = 300/10^12*1e9;    % picoJ/um^2
zeta = 0.02;            % dimensionless
epsilon = -zeta*kD;     % picoJ/um^2
% epsilon = -1;
n0 = 1;                 % fraction
kappa = 1e-19*1e12;     % picoJ
% kappa = 1e-17*1e12;     % picoJ
alpha_i = 0.01;        % fraction
zeta = epsilon*n0/kD;   % dimensionless
N = 3e3;                % number of points in quadrature
plotfigs = 0;

% param_1_vals = logspace(-21,-15, 60)*1e12;
% param_1_vals = logspace(-4,0, 40);
param_1_vals = logspace(-3,-1, 1);

for ii = 1:length(param_1_vals)
%     kappa = param_1_vals(ii);
%     zeta = param_1_vals(ii);
%     epsilon = -zeta*kD;     % picoJ/um^2
%     kD = param_1_vals(ii);
%     kappa = 3.333e-7*kD;
    sigma = param_1_vals(ii);
    d = sqrt(R^2/sigma);    % um

    % solve microplastics for initial stretch of A and B
    options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e4);
    const = [zeta, alpha_i, d, R];
    [out,fval,~,~,lam_vals,grad,hessian] = ...
        fmincon(@(y) free_phi(y, const),[alpha_i, alpha_i, 0.3], [],[],[],[],...
        [-1,-1,0],[Inf, Inf, pi/2], ...
        @(y) lipid_con_phi(y,const), options);
    alpha_A_init = out(1);
    alpha_B_init = out(2);
    phi_init = out(3);
    h_phi_init = 0;
    S_A_init = 2*pi*R^2*(1-cos(phi_init));
    S_B_init = d^2 - pi*R^2*sin(phi_init)^2;

    E_adhesion_init = epsilon*n0*S_A_init ./(1+alpha_A_init );
    E_stretch_A_init  = kD/2*(alpha_A_init .^2*S_A_init ./(1+alpha_A_init ));
    E_stretch_B_init  = kD/2*(alpha_B_init .^2*S_B_init ./(1+alpha_B_init ));
    E_total_init = E_adhesion_init+E_stretch_A_init+E_stretch_B_init;
    
    % use initial stretch to solve for minimum attachment

    const = [epsilon, n0, d, R, kD, kappa, alpha_i, N];
    opts = optimoptions(@fmincon,Algorithm="sqp");
    rng default % For reproducibility
    problem = createOptimProblem("fmincon",...
        x0 = [alpha_A_init, alpha_B_init, phi_init, h_phi_init],...
        objective = @(y) stretch_bend_min(y, const),...
        lb = [-0.1,-0.1,0,-2*d],...
        ub = [0.1, 0.1, pi/2, 2*d],...
        nonlcon = @(y) lipid_con_bend(y,const),...
        options=opts);
    % gs = GlobalSearch;
    st = alpha_i*2;
    ptmatrix = [alpha_A_init, alpha_B_init, phi_init, h_phi_init;...
        -st, -st, 5, 0;-st, st, 5, 0;st, -st, 5, 0;st, st, 5, 0;...
        -st, -st, 5, R;-st, -st, 5, -R;st, st, 5, d/2;st, st, 5, -d/2;...
        -st, -st, 85, 0;-st, st, 85, 0;st, -st, 85, 0;st, st, 85, 0;...
        -st, -st, 85, R;-st, -st, 85, -R;st, st, 85, d/2;st, st, 85, -d/2;...
        -st, -st, 45, 0;-st, st, 45, 0;st, -st, 45, 0;st, st, 45, 0;...
        -st, -st, 45, R;-st, -st, 45, -R;st, st, 45, d/2;st, st, 45, -d/2];
    tpoints = CustomStartPointSet(ptmatrix);
    rs = RandomStartPointSet('NumStartPoints',5);
    gs = MultiStart;
%     [out,fval2,exitflag,output,solutions] = run(gs,problem, {tpoints,rs});
%     [out,fval2,exitflag,output,solutions] = run(gs,problem, {rs});
    [out,fval2,exitflag,output,solutions] = run(gs,problem, {tpoints});
    
    % figure();
    % hold on
    % for jj = 1:length(solutions)
    %     plot(jj,solutions(ii).Fval, 'o')
    % end
    
    % get energies etc
    alpha_A = out(1);
    alpha_B = out(2);
    phi = out(3);
    h_phi = out(4);
    
    % get the shape of the free region
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
    r = linspace(r_phi, d/2,N);
    
    [h,C,S_B, ~, lap_h, hderiv] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi);
    
    S_A = 2*pi*R^2*(1-cos(phi));
    
    E_adhesion = epsilon*n0*S_A./(1+alpha_A);
    E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
    E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
    E_bend_B = kappa/2*2*pi*trapz(r, r.*lap_h.^2);
    E_bend_A = 4*pi*kappa*(1-cos(phi));
    
    E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;
    
    E_all(1,ii) = E;
    E_all(2,ii) = E_adhesion;
    E_all(3,ii) = E_stretch_A;
    E_all(4,ii) = E_stretch_B;
    E_all(5,ii) = E_bend_A;
    E_all(6,ii) = E_bend_B;
    
    alpha_A_vals(ii) = alpha_A;
    alpha_B_vals(ii) = alpha_B;
    phi_vals(ii) = phi;
    h_phi_vals(ii) = h_phi;
    S_A_vals(ii) = S_A;
    S_B_vals(ii) = S_B;
    Sigma_vals(ii) = Sigma;
    
    % plot of shape and nanoparticle
    if plotfigs
        figure('Position',[400,100,700,500]);
        hold on
        axis equal
        xlabel('$r$')
        ylabel('$h$')
        plot(r, h, 'displayname', 'free surface');
        t = linspace(-pi/2,pi/2,1000);
        x = cos(t)*R;
        % y = sin(t)*R+(R*cos(phi)+h(1));
        y = sin(t)*R+R*cos(phi)+h(1);
        plot(x,y, 'displayname', 'microbead')
        % plot(x, (x-sin(phi))*tan(phi)+h(1))
        % legend('location', 'se')
        annotation('textbox', [0.3,0.7,0.4,0.2], 'String',...
            [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
            sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
            sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
            sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
            sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta),...
            sprintf('$\\kappa = %0.2g k_\\mathrm{D} R^2$ \n', kappa/(kD*R^2))], ...
            'interpreter', 'latex','FontSize',16,...
            'FitBoxToText','on','LineStyle','none', 'Color','b')
        annotation('textbox', [0.6,0.7,0.4,0.2], 'String',...
            [sprintf('$\\alpha_A = %0.2g$ \n', alpha_A),...
            sprintf('$\\alpha_B = %0.2g$ \n', alpha_B),...
            sprintf('$\\phi = %0.2g^{\\circ}$ \n', rad2deg(phi)),...
            sprintf('$h_\\phi = %0.2g$ $\\mu$m \n', h_phi)], ...
            'interpreter', 'latex','FontSize',16,...
            'FitBoxToText','on','LineStyle','none', 'Color','k')
        annotation('textbox', [0.55,0.2,0.4,0.2], 'String',...
            [sprintf('$E = %0.3g$ pJ \n', E)], ...
            'interpreter', 'latex','FontSize',16,...
            'FitBoxToText','on','LineStyle','none', 'Color','r')
    end
    
end

% save('kappa_sweep_R300sig01kD300zeta02ai01.mat')

%% plotting
init_vals = param_1_vals(1:20);
on = ones(size(init_vals));

figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
ylabel('E (pJ)')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
for ii=1:6
    plot(param_1_vals, E_all(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
% plot(init_vals, on*E_total_init, 'k-.', 'linewidth', 0.5);
% plot(init_vals, on*E_adhesion_init, 'b-.', 'linewidth', 0.5);
% plot(init_vals, on*E_stretch_B_init, 'r-.', 'linewidth', 0.5);
legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
    '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
    'location', 'best')
annotation('textbox', [0.3,0.6,0.4,0.2], 'String',...
    [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
    sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
    sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
    sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
    sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta)],...
    'interpreter', 'latex','FontSize',16,...
    'FitBoxToText','on','LineStyle','none', 'Color','b')

% positive percentages
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
ylabel('Fraction of Energy')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
E_pos = E_all(1,:)-E_all(2,:);
E_pos_2 = sum(E_all(3:end,:));
E_perc = E_all./E_pos;
for ii=3:6
    plot(param_1_vals, E_perc(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
legend({'$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
    'location', 'best')

% areas
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
yyaxis left
ylabel('Area')
% plot(init_vals, on*S_B_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(init_vals, on*(S_B_init+S_A_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(param_1_vals, S_B_vals, '--', 'displayname', '$S_B$')
plot(param_1_vals, S_A_vals+S_B_vals, '-', 'displayname', '$S_\mathrm{total}$')
yyaxis right
ylabel('Area')
% plot(init_vals, on*S_A_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(param_1_vals, S_A_vals, '--', 'displayname', '$S_A$')
legend

% stretches
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
ylabel('$\alpha$')
plot(param_1_vals, alpha_B_vals, 'displayname', '$\alpha_B$')
plot(param_1_vals, alpha_A_vals, 'displayname', '$\alpha_A$')
% plot(init_vals, on*alpha_A_init, 'r-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(init_vals, on*alpha_B_init, 'b-.', 'linewidth', 0.5, 'HandleVisibility','off');
legend

% phi and h_phi
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
yyaxis left
ylabel('$\phi$')
plot(param_1_vals, rad2deg(phi_vals), 'displayname', '$\phi$')
% plot(init_vals, on*rad2deg(phi_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off')
yyaxis right
ylabel('$h_\phi$')
plot(param_1_vals, h_phi_vals, 'displayname', '$h_\phi$')
legend

%% replotting each function

figure('Position',[400,100,700,500]);
hold on
axis equal
xlabel('$r$')
ylabel('$h$')
for ii = [1,5,10,20,30,40]
    sigma = param_1_vals(ii);
    d = sqrt(R^2/sigma);    % um

    alpha_A = alpha_A_vals(ii);
    alpha_B = alpha_B_vals(ii);
    phi = phi_vals(ii);
    h_phi = h_phi_vals(ii);
    
    % get the shape of the free region
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
    r = linspace(r_phi, d/2,N);
    
    [h,C,S_B, ~, lap_h, hderiv] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi);

    h1 = plot(r, h, 'displayname', sprintf('$\\kappa = %0.2e$', kappa));
    colour = h1.Color;
    t = linspace(-pi/2,-pi/2+phi,1000);
    x = cos(t)*R;
    % y = sin(t)*R+(R*cos(phi)+h(1));
    y = sin(t)*R+R*cos(phi)+h(1);
    plot(x,y, ':', 'HandleVisibility','off', 'color', colour)
end

legend

%% extra functions
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

function f = free_phi(y, const)
    % function of free energy to minimise

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);

    alpha_A = y(1);
    alpha_B = y(2);
    phi = y(3);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = d^2 - pi*R^2*sin(phi)^2;

    f = zeta*S_A./(1+alpha_A) ...
      + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B));

end

function [c,ceq] = lipid_con_phi(y, const)

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);

    alpha_A = y(1);
    alpha_B = y(2);
    theta = y(3);

    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(theta));
    S_B = d^2 - pi*R^2*sin(theta)^2;

    c(1) = 0;

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end