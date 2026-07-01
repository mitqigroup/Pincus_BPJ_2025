c0 = 0;
lambda = -1.6;
mu = 1;

constants = [c0, lambda, mu];

p0 = 1;
initial_vals = [0,1e-5,pi/2,p0,0];

options  = odeset('Events',@surfaceEventFunc);
[s,out, se, ye, ie] = ode45(@(s,in) deriv(s,in,constants),[0,2*pi],initial_vals, options);

% [s,out] = ode45(@(s,in) deriv(s,in,constants),[0,5*pi],initial_vals);

x = out(:,1);
y = out(:,2);
Theta = out(:,3);
p = out(:,4);
b = out(:,5);

figure();
hold on
grid on
axis equal
plot(x, y);

function f = deriv(s, in, constants)
    % inputs, U = psidot
    x = in(1);
    y = in(2);
    Theta = in(3);
    p = in(4);
    b = in(5);
%     % these don't actually matter to the integration itself, just comments
%     Z = in(5);
%     A = in(6);
%     V = in(7);

    c0 = constants(1);
    lambda = constants(2);
    mu = constants(3);
    
    f(1) = cos(Theta);
    f(2) = sin(Theta);
    f(3) = p;
    f(4) = (- lambda*y.^2.*sin(Theta) - 2*p.*sin(Theta) ...
            - 2*cos(Theta).*sin(Theta)./y - b.*cos(Theta))./(2*y);
    f(5) = 2*y*lambda.*cos(Theta) + mu + p.^2 - cos(Theta).^2./y.^2 - 2*c0*p;

%     % extra constants 
%     f(5) = -sin(Theta);
%     f(6) = 2*pi*X;
%     f(7) = pi*X.^2.*sin(Theta);

    f = f'; 

end


function [ThetaEqPi, isterminal, dir] = surfaceEventFunc(s, y)
    % event function to stop integration when Theta = pi

    ThetaEqPi = y(3) + pi/2;
    isterminal = 1;
    dir = 0;

end