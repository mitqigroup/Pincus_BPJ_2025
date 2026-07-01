function [delE, alpha] = high_tension_approx(phi, constants)
% calculates the high tension approximation to the shape equations as per
% the results of Foret 2014.

epsilon = constants(1);
n0 = constants(2);
d = constants(3);
R = constants(4);
kD = constants(5);
kappa = constants(6);
alpha_i = constants(7);

% a = (1+alpha_i)*4*pi*R*sin(phi)./d^2*sqrt(kappa/kD).*(1-sqrt((1+cos(phi))/2));
% b = (1+alpha_i)*(d^2-pi*R^2*sin(phi).^2+2*pi*R^2.*(1-cos(phi)))/d^2-1;
% 
% % vals = roots([1, -(2*b+a^2), b^2])
% 
% alpha = ((2*b+a.^2)+a.*sqrt(4*b+a.^2))/2;
% or is it the other branch, alpha = ((2*b+a.^2)-a.*sqrt(4*b+a.^2))/2 ?

alpha = ones(size(phi))*alpha_i;
for ii = 1:length(phi)
    phi_val = phi(ii);
    alpha(ii) = fzero(@(alpha) alpha_diff(alpha, phi_val), [alpha_i/2, 0.1]);
end

E_bend_B = 4*pi*R*sqrt(alpha*kD*kappa).*sin(phi).*(1-sqrt((1+cos(phi))/2));
delA = 4*pi*R*sin(phi).*sqrt(kappa./(alpha*kD)).*(1-sqrt((1+cos(phi))/2));

% delA_i = 2*pi*R*sin(phi).*sqrt(alpha_i*kD/kappa).*(1-sqrt((1+cos(phi))/2));

delE = 2*pi*epsilon*n0*R^2*(1-cos(phi))./(1+alpha) ...
    + kD/2*(alpha.^2./(1+alpha).*(d^2-pi*R^2*sin(phi).^2+2*pi*R^2*(1-cos(phi))+delA)) ...
    + E_bend_B + 4*pi*kappa*(1-cos(phi)) - kD/2*d^2*(alpha_i^2/(1+alpha_i));

% can actually select lowest energy if wanted, difference is fairly minor
% alpha1 = ((2*b+a.^2)+a.*sqrt(4*b+a.^2))/2;
% alpha2 = ((2*b+a.^2)-a.*sqrt(4*b+a.^2))/2;
% % alpha = max([((2*b+a.^2)+a.*sqrt(4*b+a.^2))/2;((2*b+a.^2)-a.*sqrt(4*b+a.^2))/2]);
% % or is it the other branch, alpha = ((2*b+a.^2)-a.*sqrt(4*b+a.^2))/2 ?
% 
% E_bend_B1 = 4*pi*R*sqrt(alpha1*kD*kappa).*sin(phi).*(1-sqrt((1+cos(phi))/2));
% E_bend_B2 = 4*pi*R*sqrt(alpha2*kD*kappa).*sin(phi).*(1-sqrt((1+cos(phi))/2));
% lower_E = E_bend_B1<E_bend_B2;
% alpha = E_bend
% if 
%     alpha = alpha1;
%     E_bend_B = E_bend_B1;
% else
%     alpha = alpha2;
%     E_bend_B = E_bend_B2;
% end
% delA = 4*pi*R*sin(phi).*sqrt(kappa./(alpha*kD)).*(1-sqrt((1+cos(phi))/2));

    function diff = alpha_diff(alpha, phi)
        delA = 4*pi*R*sin(phi).*sqrt(kappa./(alpha*kD)).*(1-sqrt((1+cos(phi))/2));
        diff = alpha+1-(1+alpha_i)*(d^2-pi*R^2*sin(phi)^2+2*pi*R^2*(1-cos(phi))+delA)/d^2;
    end

end