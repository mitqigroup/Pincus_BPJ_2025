%% solve for a particular value of P, C0, sigma
R0 = 84.2;
P = 10;
C0 = 0;
Sigma = 0.88;

R = 45;
phi = pi/3;

kappa = 1e-7;
ep = 1e-3;
sig = 0.88;

psi0 = -phi;
X0 = R*sin(phi);
% psi0 = -pi/3;
% X0 = R*sin(phi);
gamma0 = 0;
Z0 = 0;
A0 = 0;
V0 = 0;

constants = [P, C0, Sigma];

% U0 = 0.2/R0;
% U0 = 0.25/R0;
% U0 = 0.245/R0;
% U0 = 6/R0;
U0 = sqrt(2*ep/kappa)+1/R;
initial_vals = [psi0,U0,X0,gamma0,Z0,A0,V0];
% initial_vals = [3.1115    1.0491    0.0316    0.0188   -0.9995   12.5647    4.1895];

options  = odeset('Events',@surfaceEventFunc);

[s,out, se, ye, ie] = ode45(@(s,y) deriv(s,y,constants),[0,pi*5],initial_vals, options);

psi = out(:,1);
U = out(:,2);
X = out(:,3);
gamma = out(:,4);
Z = out(:,5);
A = out(:,6);
V = out(:,7);

final_vals = [psi(end), U(end), X(end), gamma(end), Z(end), A(end), V(end)]

% %%
figure();
hold on
grid on
axis equal
plot(X, Z, 'r-', 'displayname', 'Variational');
t = linspace(-pi,pi,100);
x = R*sin(t);
z = R*cos(t)+R*cos(phi);
psi_test = t;
plot(x,z, 'r--', 'displayname', 'Sphere')

%% Fig 16 Lipowsky 1991
% this is using many values of U0 and determining X(S_end) for each value
% of U0. Now we want to update with a given initial curvature and X

% psi = y(1);
% U = y(2);
% X = y(3);
% gamma = y(4);
% % these don't actually matter to the integration itself, just comments
% Z = y(5);
% A = y(6);
% V = y(7);

R = 0.3;
% phi = pi/1000;
phi = 2*pi/3;

psi0 = -phi;
X0 = R*sin(phi)+1e-10;
gamma0 = 0;
Z0 = 0;
A0 = 0;
V0 = 0;

R0 = 1;
P = 1;
C0 = 0;
Sigma = -1.1*P^(2/3);
% Sigma = -0.5*P*R0;
% U0_vals = (-2:0.01:3)/P^(-1/3);
U0_vals = linspace(-10,10,500);
maxN = 3;
X_end = nan(length(U0_vals), maxN);
S_end = nan(size(X_end));
U_end = nan(size(X_end));

% P = 1;
% C0 = 0;
% Sigma = -1.1*P^(2/3);
% U0_vals = (-2:0.01:3)/P^(-1/3);
% X_end = nan(size(U0_vals));
% S_end = nan(size(U0_vals));
% U_end = nan(size(U0_vals));

constants = [P, C0, Sigma];

for ii = 1:length(U0_vals)
    U0 = U0_vals(ii);

    initial_vals = [psi0,U0,X0,gamma0,Z0,A0,V0];
    options  = odeset('Events',@surfaceEventFunc);
    
    [s,out, se, ye, ie] = ode23(@(s,y) deriv(s,y,constants),[0,5*pi],initial_vals, options);
    
    size_ye = size(ye);
    for nn = 1:min(maxN, size_ye(1))
        X_end(ii,nn) = ye(nn,3);
        S_end(ii,nn) = se(nn);
        U_end(ii,nn) = ye(nn,2);
        Z_end(ii,nn) = ye(nn,5);
    end
end

% U0 vs X
f1 = figure('Position', [500 100 700 500]);
hold on
axes1 = gca;
fsize=20;
% axes1.YScale = 'log';
% axes1.XScale = 'log';
hold(axes1,'on');
box(axes1,'on');
grid(axes1, 'on');
set(axes1,'FontSize',16,'LineWidth',2,'TickLength',[0.015 0.025]);
xlabel('$U(0) \bar{P}^{-1/3}$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$X(S^{(n)}_1) \bar{P}^{1/3}$', 'Interpreter', 'latex', 'FontSize', 20)
% title(strcat('$N = $ ', num2str(N), ', $h^* = $',hstarstring),...
%     'Interpreter', 'latex', 'FontSize', 24)

symbols = ["x", "o", "^"];
colours = ["#27E0D8", "#648FFF", "#1D0C6B", "#DC267F", "#FE6100", ...
    "#FFB000", "#016D24"];
for nn = 1:maxN
    X_plot = X_end(:,nn);
%     X_plot = X_end;
%     X_plot(X_plot==0) = nan;
    plot(U0_vals*P^(-1/3), X_plot*P^(1/3), 'kx', 'color', colours(nn), 'Marker', symbols(nn),...
        'displayname', sprintf('$n = %d$', nn))
end

legend('location','best', 'interpreter', 'latex')

%% run forwards/backwards with function

[selection,~] = ginput(1);
index = find(U0_vals*P^(-1/3)<selection, 1, 'last');

U0_star = selection;
U1_star = U_end(index,1);
S_star = S_end(index,1);
Z_star = Z_end(index,1);

inp0 = [U0_star, U1_star, S_star, Z_star];

options = optimset('TolFun', 1e-10, 'TolX', 1e-10, 'MaxFunEvals', 1e5);

vals = fminsearch(@(inp) shoot_forwards_back(inp, constants, R, phi), inp0, ...
    options);

% vals = fmincon(@(inp) shoot_forwards_back(inp, constants, R, phi), inp0, [], []);

U0_tilde = vals(1);
U1_tilde = vals(2);
S_tilde = vals(3);
Z_tilde = vals(4);

% forwards integration
s_bar = S_tilde/2;
initial_vals = [psi0,U0_tilde,X0,gamma0,Z0,A0,V0];
[s,out] = ode78(@(s,y) deriv(s,y,constants),[0,s_bar*2],initial_vals);

psi_f = out(:,1);
U_f = out(:,2);
X_f = out(:,3);
gamma_f = out(:,4);
Z_f = out(:,5);
A_f = out(:,6);
V_f = out(:,7);

figure();
hold on
axis equal
plot(X_f, Z_f, 'b-');
t = linspace(-pi,pi,100);
x = R*sin(t);
z = R*cos(t)+R*cos(phi);
psi_test = t;
plot(x,z, 'r--', 'displayname', 'Sphere')

% backwards integration
initial_vals = [pi,U1_tilde,1e-10,0,Z_tilde,0,0];
[s,out] = ode78(@(s,y) deriv(s,y,constants),[S_tilde,0],initial_vals);

psi_b = out(:,1);
U_b = out(:,2);
X_b = out(:,3);
gamma_b = out(:,4);
Z_b = out(:,5);
A_b = out(:,6);
V_b = out(:,7);

plot(X_b, Z_b, 'r-');

total_vol = V_f(end) - V_b(end);
total_area = A_f(end) - A_b(end);

reduced_volume = 6*sqrt(pi)*total_vol/total_area^(3/2);

% %% run forwards/backwards with those values
% 
% [selection,~] = ginput(1);
% index = find(U0_vals*P^(-1/3)<selection, 1, 'last');
% 
% U0_star = selection;
% U1_star = U_end(index,1);
% S_star = S_end(index,1);
% Z_star = Z_end(index,1);
% 
% % forwards integration
% s_bar = S_star/2;
% initial_vals = [0,U0_star,1e-10,0,0,0,0];
% [s,out] = ode78(@(s,y) deriv(s,y,constants),[0,s_bar],initial_vals);
% 
% psi_bar = out(end,1);
% U_bar = out(end,2);
% X_bar = out(end,3)
% gamma_bar = out(end,4);
% 
% psi = out(:,1);
% U = out(:,2);
% X = out(:,3);
% gamma = out(:,4);
% Z = out(:,5);
% A = out(:,6);
% V = out(:,7);
% 
% figure();
% hold on
% axis equal
% plot(X, Z, 'b-');
% 
% % backwards integration
% initial_vals = [pi,U1_star,1e-10,0,Z_star,0,0];
% [s,out] = ode78(@(s,y) deriv(s,y,constants),[S_star,s_bar],initial_vals);
% 
% psi_hat = out(end,1);
% U_hat = out(end,2);
% X_hat = out(end,3)
% gamma_hat = out(end,4);
% 
% psi = out(:,1);
% U = out(:,2);
% X = out(:,3);
% gamma = out(:,4);
% Z = out(:,5);
% A = out(:,6);
% V = out(:,7);
% 
% plot(X, Z, 'rx');
% 
% 
% 
% %% forwards/backwards testing
% R0 = 1;
% P = 1;
% C0 = 0;
% Sigma = -1/2*P*R0;
% 
% constants = [P, C0, Sigma];
% 
% U0_star = 2/R0;
% s_star = pi;
% s_bar = 3;
% 
% options  = odeset('RelTol',1e-8);
% 
% % forward integration
% initial_vals = [0,U0_star,1e-2,0,0,0,0];
% [s,out] = ode78(@(s,y) deriv(s,y,constants),[0,s_bar],initial_vals, options);
% 
% final_vals = out(end,:)
% 
% psi = out(:,1);
% U = out(:,2);
% X = out(:,3);
% gamma = out(:,4);
% Z = out(:,5);
% A = out(:,6);
% V = out(:,7);
% 
% figure();
% hold on
% axis equal
% plot(X, Z, 'b-');
% 
% % %backwards integration
% [s,out] = ode78(@(s,y) deriv(s,y,constants),[s(end:-1:1)],final_vals, options);
% 
% % should be the same as initial conditions to within error
% out(end,:)
% 
% psi = out(:,1);
% U = out(:,2);
% X = out(:,3);
% gamma = out(:,4);
% Z = out(:,5);
% A = out(:,6);
% V = out(:,7);
% 
% % figure();
% hold on
% axis equal
% plot(X, Z, 'ro');
% 
% %% forwards/backwards real
% R0 = 1;
% P = 1;
% C0 = 0;
% Sigma = -1/2*P*R0;
% 
% constants = [P, C0, Sigma];
% 
% U0_star = 2/R0;
% s_star = pi;
% s_bar = 1;
% s_hat = s_star-s_bar;
% U1_star = 2/R0;
% 
% options  = odeset('RelTol',1e-8);
% 
% % forward integration
% initial_vals = [0,U0_star,1e-10,0,2,0,0];
% [s,out] = ode78(@(s,y) deriv(s,y,constants),[0,s_bar],initial_vals, options);
% 
% psi_bar = out(end,1);
% U_bar = out(end,2);
% X_bar = out(end,3)
% gamma_bar = out(end,4);
% 
% psi = out(:,1);
% U = out(:,2);
% X = out(:,3);
% gamma = out(:,4);
% Z = out(:,5);
% A = out(:,6);
% V = out(:,7);
% 
% figure();
% hold on
% axis equal
% plot(X, Z, 'b-');
% 
% % %backwards integration
% initial_vals = [pi,U1_star,1e-10,0,0,0,0];
% [s,out] = ode78(@(s,y) deriv(s,y,constants),[s_star, s_bar],initial_vals, options);
% 
% psi_hat = out(end,1);
% U_hat = out(end,2);
% X_hat = out(end,3)
% gamma_hat = out(end,4);
% 
% psi = out(:,1);
% U = out(:,2);
% X = out(:,3);
% gamma = out(:,4);
% Z = out(:,5);
% A = out(:,6);
% V = out(:,7);
% 
% % figure();
% % hold on
% % axis equal
% plot(X, Z, 'r-');
% 
% % figure();
% % hold on
% % axes1 = gca;
% % % axes1.YScale = 'log';
% % % axes1.XScale = 'log';
% % plot(X, -U, 'r-');

%% functions

function f = deriv(s, y, constants)
    % inputs, U = psidot
    psi = y(1);
    U = y(2);
    X = y(3);
    gamma = y(4);
    % these don't actually matter to the integration itself, just comments
    Z = y(5);
    A = y(6);
    V = y(7);

    P = constants(1);
    C0 = constants(2);
    Sigma = constants(3);
    
    f(1) = U;
    f(2) = -U./X.*cos(psi) + cos(psi).*sin(psi)./X.^2 + ...
           gamma./X.*sin(psi) + P*X/2.*cos(psi);
%     f(2) = -U./X.*cos(psi) + cos(psi).*sin(psi)./X.^2 + ...
%            gamma./(4*X).*sin(psi) + P*X/8.*cos(psi);
%     f(2) = -U./X.*cos(psi) + cos(psi).*sin(psi)./X.^2 + ...
%            gamma./(2*X).*sin(psi) + P*X/4.*cos(psi);
    f(3) = cos(psi);
    f(4) = (C0-U).^2/2 - sin(psi).^2./(2*X.^2) + P*X.*sin(psi) + Sigma;
%     f(4) = (C0-2*U).^2/2 - 2*sin(psi).^2./(X.^2) + P*X.*sin(psi) + Sigma;
%     f(4) = (C0-U).^2 - sin(psi).^2./(X.^2) + P*X.*sin(psi) + Sigma;

    % extra constants 
    f(5) = -sin(psi);
    f(6) = 2*pi*X;
    f(7) = pi*X.^2.*sin(psi);

    f = f'; 

end


function [psiEqPi, isterminal, dir] = surfaceEventFunc(s, y)
    % event function to stop integration when psi = pi

    psiEqPi = y(1) - pi;
    isterminal = 0;
    dir = 0;

end

function diffs = shoot_forwards_back(inputs, constants, R, phi)
    % we're going to minimise the difference between the output ends

%     R = 0.3;
%     phi = pi/6;
    
    psi0 = -phi;
    X0 = R*sin(phi)+1e-10;

    U0_star = inputs(1);
    U1_star = inputs(2);
    S_star = inputs(3);
    Z_star = inputs(4);
    
    % forwards integration
    s_bar = S_star/2;
    initial_vals = [psi0,U0_star,X0,0,0,0,0];
    [s,out_forward] = ode78(@(s,y) deriv(s,y,constants),[0,s_bar],initial_vals);
    
%     psi_bar = out(end,1);
%     U_bar = out(end,2);
%     X_bar = out(end,3);
%     gamma_bar = out(end,4);
%     Z_bar = out(end,5);
    
    % backwards integration
    initial_vals = [pi,U1_star,1e-10,0,Z_star,0,0];
    [s,out_backward] = ode78(@(s,y) deriv(s,y,constants),[S_star,s_bar],initial_vals);
    
%     psi_hat = out(end,1);
%     U_hat = out(end,2);
%     X_hat = out(end,3);
%     gamma_hat = out(end,4);
%     Z_hat = out(end,5);

    diffs = sum((out_forward(end,1:5)-out_backward(end,1:5)).^2);
%     diffs = sum((out_forward(end,1:3)-out_backward(end,1:3)).^2);

%     diffs = (psi_bar-psi_hat).^2 + ...
%         (U_bar-U_hat).^2 + ...
%         (X_bar-X_hat).^2 + ...
%         (gamma_bar-gamma_hat).^2;

end