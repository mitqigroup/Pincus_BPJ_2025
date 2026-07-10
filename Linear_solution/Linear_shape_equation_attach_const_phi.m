% clear variables

% constants
R = 0.5;
sigma = 3.5e-3;
d = sqrt(R^2/sigma);    % um
% d = 20;    % um
phi = pi/12;
kD  = 300/10^12*1e9;    % picoJ/um^2
zeta = 2.7e-6;            % dimensionless
epsilon = -zeta*kD;     % picoJ/um^2
n0 = 1;                 % fraction
kappa = 1e-19*1e12;     % picoJ
kappa_bar = kappa/(kD*R^2);
% gamma = 0.5*R;            % um
% epsilon = -kappa/gamma^2;
% kappa = 1e-17*1e12;     % picoJ
% zeta = -epsilon*n0/kD;   % dimensionless
% alpha_i = sqrt(-zeta/2);        % fraction
alpha_i = 0.01;
% alpha_i = 0;
% ten_param=-d*epsilon*n0/sqrt(kD*kappa)
N = 3e3;                % number of points in quadrature

savedata = 0;
plotfigs = 0;
filename = '300nm_ai1kap19zeta002.mat';

phi_vals = flip(deg2rad(linspace(0.001,0.1,50)));

% phi_vals = deg2rad(30);

E_all = zeros(6,length(phi_vals));

tic
for ii = 1:length(phi_vals)
phi = phi_vals(ii);
rad2deg(phi)

% if ii==1
% %     % solve for initial stretch of A and B
% %     options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e4);
% %     const = [zeta, alpha_i, d, R, phi];
% %     [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
% %         fmincon(@(y) free_phi(y, const),[alpha_i, alpha_i], [],[],[],[],...
% %         [-1,-1],[Inf, Inf], ...
% %         @(y) lipid_con_phi(y,const), options);
% %     
% %     alpha_A_init = out(1)
% %     alpha_B_init = out(2)
% %     h_phi_init = 0;
% %     
% %     lambda_stretch = lam_vals.eqnonlin;
% %     
% %     alpha_A_test = sqrt(1+2*(lambda_stretch+zeta))-1;
% %     alpha_B_test = sqrt(1+2*lambda_stretch)-1;
% 
% % 
%     
% %     alpha_B_init = alpha_i+0.001;
% %     alpha_A_init = 0.05;
% %     h_phi_init = 0;
% else
%     alpha_A_init = alpha_A;
%     alpha_B_init = alpha_B;
%     h_phi_init = h_phi;
% end

if ii==1
    alpha_A_init = (alpha_i-zeta)/(1+zeta);
    alpha_B_init = alpha_i;
    h_phi_init = 0;
else
    alpha_A_init = alpha_A;
    alpha_B_init = alpha_B;
    h_phi_init = h_phi;
end

% alpha_A_init = -0.02;
% alpha_B_init = -0.01;
% h_phi_init = 0.3;

%% use initial stretch to solve for minimum attachment
% options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e3);
% options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp');
% const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi];
% [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
%     fmincon(@(y) stretch_bend_min(y, const),[alpha_A_init, alpha_B_init, h_phi_init],...
%     [],[],[],[],[-0.1,-0.1,-Inf],[0.1, 0.1, Inf], ...
%     @(y) lipid_con_bend(y,const), options);

const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi];
[out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_linear_minimum(const, [alpha_A_init, alpha_B_init, h_phi_init]);

alpha_A = out(1);
alpha_B = out(2);
h_phi = out(3);

lambda_stretch = lam_vals.eqnonlin;
% 
% [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
%     fmincon(@(y) stretch_bend_min_explicit_alpha_A(y, const),[alpha_B_init, h_phi_init],...
%     [],[],[],[],[-0.1,-Inf],[0.1, Inf],[],options);

% options = optimset('Display','iter','PlotFcns',@optimplotfval);
% out = fminsearch(@(y) stretch_bend_min_explicit_alpha_A(y, const),[0, 0], options);
% 
% alpha_B = out(1)
% h_phi = out(2);

% alpha_A_test = sqrt(1+2*(lambda_stretch+zeta))-1;
% alpha_B_test = sqrt(1+2*lambda_stretch)-1;

% get the shape of the free region
Sigma = kD*alpha_B;
lambda = sqrt(kappa/Sigma);
r_phi = sin(phi)*R;
r = linspace(r_phi, d/2,N);

[h,C,S_B, ~, hderiv, lap_h] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi);

S_A = 2*pi*R^2*(1-cos(phi));
alpha_A_test = S_A*(d^2/(1+alpha_i)-S_B/(1+alpha_B))^(-1)-1;
alpha_A_vals_test(ii) = alpha_A_test;

S_B_test(ii) = -pi/2*r_phi*tan(phi)^2*(r_phi...
    -(besselk(0,r_phi/lambda,1)*(2*lambda*besselk(1,r_phi/lambda,1)...
    +r_phi*besselk(0,r_phi/lambda,1))/(besselk(1,r_phi/lambda,1)^2)))...
    +(d^2-pi*r_phi^2);
% alpha_A = S_A*(d^2/(1+alpha_i)-S_B/(1+alpha_B))^(-1)-1

E_adhesion = epsilon*n0*S_A./(1+alpha_A);
E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
E_bend_B = kappa/2*2*pi*trapz(r, r.*lap_h.^2);
E_bend_test = pi/2*kD*alpha_B.*R^2.*tan(phi).^2.*sin(phi).^2.*...
    (besselk(1,R./lambda.*sin(phi),1).^2-...
     besselk(0,R./lambda.*sin(phi),1).^2)./...
     besselk(1,R./lambda.*sin(phi),1).^2;
E_bend_A = 4*pi*kappa*(1-cos(phi));

E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;

initial_E = alpha_i^2*kD*d^2/(1+alpha_i)/2;

E_all(1,ii) = E-initial_E;
E_all(2,ii) = E_adhesion;
E_all(3,ii) = E_stretch_A;
E_all(4,ii) = E_stretch_B-initial_E;
E_all(5,ii) = E_bend_A;
E_all(6,ii) = E_bend_B;

alpha_A_vals(ii) = alpha_A;
alpha_B_vals(ii) = alpha_B;
h_phi_vals(ii) = h_phi;
S_A_vals(ii) = S_A;
S_B_vals(ii) = S_B;
Sigma_vals(ii) = Sigma;

% hderiv = C(1)./r+C(3)/lambda*besseli(1,r/lambda)-C(4)/lambda*besselk(1,r/lambda);
% lap_h = C(3)/lambda^2*besselk(0,r/lambda)+C(4)/lambda^2*besselk(0,r/lambda);

%% plot of shape and nanoparticle
% subplot(1,2,1);
if plotfigs
figure('Position',[400,100,800,600]);
hold on
% axis equal
xlabel('$r$')
ylabel('$h$')
plot(r, h-h(1), 'displayname', 'free surface');
t = linspace(-pi/2,pi/2,1000);
x = cos(t)*R;
% y = sin(t)*R+(R*cos(phi)+h(1));
y = sin(t)*R+R*cos(phi);
plot(x,y, 'displayname', 'microbead')
% plot(x, (x-sin(phi))*tan(phi)+h(1))
% legend('location', 'se')
annotation('textbox', [0.3,0.7,0.4,0.2], 'String',...
    [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
    sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
    sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
    sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
    sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta),...
    sprintf('$\\kappa = %0.2g$ pJ \n', kappa)], ...
    'interpreter', 'latex','FontSize',16,...
    'FitBoxToText','on','LineStyle','none', 'Color','b')
annotation('textbox', [0.55,0.7,0.4,0.2], 'String',...
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
% % annotation('textbox', [0.55,0.2,0.4,0.2], 'String',...
% %     [sprintf('$E = %0.2g$ \n', E),...
% %     sprintf('$E_\\mathrm{adhesion} = %0.1f\\%%$ \n', E_adhesion/E*100), ...
% %     sprintf('$E_\\mathrm{stretch,A} = %0.1f\\%%$ \n', E_stretch_A/E*100), ...
% %     sprintf('$E_\\mathrm{stretch,B} = %0.1f\\%%$ \n', E_stretch_B/E*100), ...
% %     sprintf('$E_\\mathrm{bend,A} = %0.1f\\%%$ \n', E_bend_A/E*100), ...
% %     sprintf('$E_\\mathrm{bend,B} = %0.1f\\%%$ \n', E_bend_B/E*100)], ...
% %     'interpreter', 'latex','FontSize',16,...
% %     'FitBoxToText','on','LineStyle','none', 'Color','g')
% 
% % subplot(1,2,2);
% axes('Position',[.5 .3 .4 .4])
% % figure();
% hold on
% ylabel('Energy (pJ)')
% X = categorical({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
%     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
%     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'});
% X = reordercats(X,{'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
%     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
%     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'});
% bar(X, E_all, 'barwidth', 1)

% set(gca,'xticklabel',{'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
%     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
%     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'});
% bar([-0.5:4.5],E_all, 'barwidth', 1)

% plot of derivative
% figure();
% hold on
% plot(r,hderiv)
% plot(r, tan(phi)*ones(size(r)))
end

end
toc

if savedata
    save(filename);
end

%%
% % close all
% if ~ishandle(1)
%     f1 = figure('Position',[400,100,700,500]);
% else
%     figure(f1);
% end
hold on
% xlim([0,90])
xlabel('$\phi$')
ylabel('$\Delta E$')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
for ii=1:1
%     plot(rad2deg(phi_vals), E_all(ii,:), ...
%         strcat(colours(ii),lines(ii)))
    plot(rad2deg(phi_vals), E_all(ii,:), '-');
end
% annotation('textbox', [0.5,0.7,0.4,0.2], 'String',...
%     [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
%     sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
%     sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
%     sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
%     sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta),...
%     sprintf('$\\kappa = %0.2g k_\\mathrm{D} R^2$ \n', kappa/(kD*R^2))], ...
%     'interpreter', 'latex','FontSize',16,...
%     'FitBoxToText','on','LineStyle','none', 'Color','b')
% legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
%     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
%     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
%     'location', 'nw')

% 
% %slopes
% if ~ishandle(2)
%     f2 = figure('Position',[400,100,700,500]);
% else
%     figure(f2);
% end
% hold on
% % xlim([0,90])
% xlabel('$\phi$')
% ylabel('$\partial E/\partial \phi$')
% lines = ["-", ":", ":", "--", ":", "--"];
% colours = ['k', 'b', 'r', 'r', 'g', 'g'];
% for ii=1:6
%     plot(rad2deg((phi_vals(1:end-1)+phi_vals(2:end))/2),...
%         diff(E_all(ii,:))./diff(rad2deg(phi_vals)), ...
%         strcat(colours(ii),lines(ii)))
% end
% % annotation('textbox', [0.5,0.7,0.4,0.2], 'String',...
% %     [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
% %     sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
% %     sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
% %     sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
% %     sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta),...
% %     sprintf('$\\kappa = %0.2g k_\\mathrm{D} R^2$ \n', kappa/(kD*R^2))], ...
% %     'interpreter', 'latex','FontSize',16,...
% %     'FitBoxToText','on','LineStyle','none', 'Color','b')
% % legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
% %     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
% %     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
% %     'location', 'nw')
% 
% 
% % % areas
% % figure();
% % hold on
% % xlim([0,90]);
% % xlabel('$\phi$')
% % yyaxis left
% % ylabel('Area')
% % plot(rad2deg(phi_vals), S_B_vals, '--', 'displayname', '$S_B$')
% % plot(rad2deg(phi_vals), S_A_vals+S_B_vals, '-', 'displayname', '$S_\mathrm{total}$')
% % yyaxis right
% % ylabel('Area')
% % plot(rad2deg(phi_vals), S_A_vals, '--', 'displayname', '$S_A$')
% % legend
% % 
% % % stretches
% figure();
% hold on
% % xlim([0,90]);
% xlabel('$\phi$')
% ylabel('$\alpha$')
% plot(rad2deg(phi_vals), alpha_B_vals, 'displayname', '$\alpha_B$')
% plot(rad2deg(phi_vals), alpha_A_vals, 'displayname', '$\alpha_A$')
% % plot(rad2deg(phi_vals), alpha_A_vals_test, 'displayname', '$\alpha_A test$')
% legend

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
    phi = const(9);
    
    alpha_A = y(1);
    alpha_B = y(2);
    h_phi = y(3);
    
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

function f = stretch_bend_min_explicit_alpha_A(y, const)
    % function of free energy to minimise

    epsilon = const(1);
    n0 = const(2);
    d = const(3);
    R = const(4);
    kD = const(5);
    kappa = const(6);
    alpha_i = const(7);
    N = const(8);
    phi = const(9);
    
    alpha_B = y(1);
    h_phi = y(2);
    
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
    r = linspace(r_phi, d/2,N);
    
    [~,~,S_B, ~, lap_h, ~] = free_shape_linear_fixed_h(...
        r, r_phi, d, phi, kappa, Sigma, h_phi);

    S_A = 2*pi*R^2*(1-cos(phi));

    alpha_A = S_A*(d^2/(1+alpha_i)-S_B/(1+alpha_B))^(-1)-1;
    
    % stretching, adhesion and bending energy
    f = epsilon*n0*S_A./(1+alpha_A) ...
      + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
      + kappa/2*2*pi*trapz(r, r.*lap_h.^2) + 4*pi*kappa*(1-cos(phi));

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
    phi = const(9);
    
    alpha_A = y(1);
    alpha_B = y(2);
    h_phi = y(3);
    
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
    phi = const(5);

    alpha_A = y(1);
    alpha_B = y(2);

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
    phi = const(5);

    alpha_A = y(1);
    alpha_B = y(2);

    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = d^2 - pi*R^2*sin(phi)^2;

    c(1) = 0;

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end

function [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_linear_minimum(constants, inputs)

    epsilon = constants(1);
    n0 = constants(2);
    d = constants(3);
    R = constants(4);
    kD = constants(5);
    kappa = constants(6);
    alpha_i = constants(7);
    N = constants(8);
    phi = constants(9);
    
    alpha_A_init = inputs(1);
    alpha_B_init = inputs(2);
    h_phi_init = inputs(3);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = 0;
    S_i = d^2;

    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;
    h_phi = h_phi_init;

%     options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp');

    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
        'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12);

    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init, h_phi_init],...
        [],[],[],[],[-0.1,-0.1,-Inf],[0.1, 0.1, Inf], ...
        @constraint, options);

    function f = objective(y_obj)

        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        h_phi = y_obj(3);
        
        Sigma = kD*alpha_B;
        lambda = sqrt(kappa/Sigma);
        r_phi = sin(phi)*R;
        r = linspace(r_phi, d/2,N);
        
        [~,~,S_B, ~,  ~, lap_h] = free_shape_linear_fixed_h(...
            r, r_phi, d, phi, kappa, Sigma, h_phi);

        % stretching, adhesion and bending energy
        f = epsilon*n0*S_A./(1+alpha_A) ...
          + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
          + kappa/2*2*pi*trapz(r, r.*lap_h.^2) + 4*pi*kappa*(1-cos(phi));

    end

    function [c,ceq]= constraint(~)
    
        c(1) = 0;
        ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
    
    end

end