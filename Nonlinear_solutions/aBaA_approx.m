function [delE, alpha] = aBaA_approx(phi, constants)
% calculates the high tension approximation to the shape equations as per
% the results of Foret 2014.

epsilon = constants(1);
n0 = constants(2);
d = constants(3);
R = constants(4);
kD = constants(5);
kappa = constants(6);
alpha_i = constants(7);

alpha = (1+alpha_i)*(d^2-pi*R^2*sin(phi).^2+2*pi*R^2.*(1-cos(phi)))/d^2-1;

E_bend_B = 0;
delA = 0;

delE = 2*pi*epsilon*n0*R^2*(1-cos(phi))./(1+alpha) ...
    + kD/2*(alpha.^2./(1+alpha).*(d^2-pi*R^2*sin(phi).^2+2*pi*R^2*(1-cos(phi))+delA)) ...
    + E_bend_B + 4*pi*kappa*(1-cos(phi)) - kD/2*d^2*(alpha_i^2/(1+alpha_i));

end