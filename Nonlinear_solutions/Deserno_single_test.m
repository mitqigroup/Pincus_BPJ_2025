% constants and inputs
R = 1;
phi = deg2rad(30);
Sig_bar = 1;
lambda = R/sqrt(Sig_bar);
const = [lambda];
% psi_dot_init = -1.275977787;
psi_dot_init = -1.27598;
% psi_dot_init = -1.27597;
% psi_dot_init = -2;
% psi_dot_init = psi_dot_init*1.0001;
% psi_dot_init = psi_dot_init*0.99999999;

% intial conditions and boundary conditions
r_phi = R*sin(phi);
psi_init = phi;
p_psi_init = 2*r_phi*(psi_dot_init+sin(psi_init)/r_phi);
p_r_init = -(p_psi_init^2/(4*r_phi)-p_psi_init*sin(psi_init)/r_phi...
    -2*r_phi/lambda^2)/cos(psi_init);

% solving ODE
init_vals = [psi_init, r_phi, R*sin(phi), p_psi_init, p_r_init, 0,0,0,0];
event_func = @(s,y) myEventFcn(s,y,const);
%     options = odeset('Events', @myEventFcn);
options = odeset('Events', event_func, 'RelTol',1e-12,'AbsTol',1e-13);
% options = odeset();
out = ode45(@(s, y) hamilton(s, y,const),[0,30],init_vals,...
    options);

% plot solution
% r = out.y(2,:);
% h = out.y(3,:);
yvals = deval(out, linspace(0,30,1000));
r = yvals(2,:);
h = yvals(3,:);
% figure('Position',[400,100,800,600]);
hold on
axis equal
xlabel('$r$')
ylabel('$h$')
plot(r, h-h(1), 'b-', 'displayname', 'free surface');
t = linspace(0,pi,1000);
x = sin(t)*R;
% y = sin(t)*R+(R*cos(phi)+h(1));
y = cos(t)*R+R*cos(phi);
% plot(x,y,'r:', 'displayname', 'microbead')

function derivs = hamilton(s, y, const)
    psi = y(1);
    r = y(2);
    h = y(3);
    p_psi = y(4);
    p_r = y(5);
    p_h = y(6);
    E_bend = y(7);
    A = y(8);
    E_free = y(9);

    lambda = const(1);

    f(1) = p_psi/(2*r)-sin(psi)/r;
    f(2) = cos(psi);
    f(3) = sin(psi);
    f(4) = cos(psi)*(p_psi/r-p_h)...
        +sin(psi)*p_r;
    f(5) = p_psi/r*(p_psi/(4*r)-sin(psi)/r)...
        +2/lambda^2;
    f(6) = 0;
    f(7) = p_psi^2/(4*r);
    f(8) = r;
    f(9) = f(7) + r*2/lambda^2*(1-cos(psi));

    derivs = f';
end

function [value,isterminal,direction] = myEventFcn(s,y,const)
    psi = y(1);
    r = y(2);
    h = y(3);
    p_psi = y(4);
    p_r = y(5);
    p_h = y(6);
    E_bend = y(7);
    A = y(8);
    E_free = y(9);

    lambda = const(1);
    
    value = [psi, psi-pi];
%     isterminal = [1,1];
    isterminal = [0,0];
    direction = [0,0];
end