phi_vals = linspace(0.001, pi/4, 50);

for ii = 1:length(phi_vals)

phi = phi_vals(ii);
% constants
R = 0.1;                  % um
sigma = 0.011;          % surface fraction
d = sqrt(R^2/sigma);    % um
kD  = 300/10^12*1e9;    % picoJ/um^2
zeta = 0.02;            % dimensionless
epsilon = -zeta*kD;     % picoJ/um^2
% epsilon = -1;
n0 = 1;                 % fraction
kappa = 1e-19*1e12;     % picoJ
alpha_i = 0.013;        % fraction
zeta = epsilon*n0/kD;   % dimensionless
N = 3e3;                % number of points in quadrature

% solve for initial stretch of A and B
options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e4);
const = [zeta, alpha_i, d, R, phi];
[out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    fmincon(@(y) free_phi(y, const),[0.01, 0.01], [],[],[],[],...
    [-1,-1],[Inf, Inf], ...
    @(y) lipid_con_phi(y,const), options);

alpha_A_init = out(1);
alpha_B_init = out(2);

lambda_stretch = lam_vals.eqnonlin;

alpha_A_test = sqrt(1+2*(lambda_stretch+zeta))-1;
alpha_B_test = sqrt(1+2*lambda_stretch)-1;

%% use initial stretch to solve for minimum attachment
options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e3);
const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi];
[out,Energy(ii),exitflag,output,lam_vals,grad,hessian] = ...
    fmincon(@(y) stretch_bend_min(y, const),[alpha_A_init, alpha_B_init], [],[],[],[],...
    [-1,-1],[Inf, Inf], ...
    @(y) lipid_con_bend(y,const), options);

alpha_A = out(1);
alpha_B = out(2);

lambda_stretch = lam_vals.eqnonlin;

% alpha_A_test = sqrt(1+2*(lambda_stretch+zeta))-1;
% alpha_B_test = sqrt(1+2*lambda_stretch)-1;

% get the shape of the free region
Sigma = kD*alpha_B;
lambda = sqrt(kappa/Sigma);
r_phi = sin(phi)*R;
r = linspace(r_phi, d/2,N);

[h,C,A] = free_shape_linear(r, R, d, phi, kappa, Sigma, -epsilon*n0/(1+alpha_A));

hderiv = C(1)./r+C(3)/lambda*besseli(1,r/lambda)-C(4)/lambda*besselk(1,r/lambda);
lap_h = C(3)/lambda^2*besselk(0,r/lambda)+C(4)/lambda^2*besselk(0,r/lambda);

% plot of shape and nanoparticle
% figure();
% hold on
% axis equal
% plot(r, h);
% t = linspace(-pi/2,pi/2,1000);
% x = cos(t)*R;
% % y = sin(t)*R+(R*cos(phi)+h(1));
% y = sin(t)*R+R*cos(phi)+h(1);
% plot(x,y)
% plot(x, (x-sin(phi))*tan(phi)+h(1))

% plot of derivative
% figure();
% hold on
% plot(r,hderiv)
% plot(r, tan(phi)*ones(size(r)))

end

figure()
plot(phi_vals, Energy)

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
    
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
    r = linspace(r_phi, d/2,N);
    
    [h,C,S_B] = free_shape_linear(r, R, d, phi, kappa, Sigma, -epsilon*n0/(1+alpha_A));
    
%     hderiv = C(1)./r+C(3)/lambda*besseli(1,r/lambda)-C(4)/lambda*besselk(1,r/lambda);
    lap_h = C(3)/lambda^2*besselk(0,r/lambda)+C(4)/lambda^2*besselk(0,r/lambda);
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
    phi = const(9);
    
    alpha_A = y(1);
    alpha_B = y(2);
    
    Sigma = kD*alpha_B;
    r_phi = sin(phi)*R;
    
    r = linspace(r_phi, d/2,N);
    [h,C,S_B] = free_shape_linear(r, R, d, phi, kappa, Sigma, -epsilon*n0/(1+alpha_A));
    
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