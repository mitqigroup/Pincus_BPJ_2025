function [h, C, A] = free_shape_linear(r, R, d, phi, kappa, Sigma, w)
% solves the shape equation for the free surface of a section of membrane
% bound at one end to a sphere of radius R at the point phi radians from
% the bottom of the sphere, with a flat surface at d/2. The shape is
% characterised by the length lambda, which is equal to the square root of
% the membrane bending energy kappa divided by the tension Sigma. Also
% outputs the total area of the membrane, as well as the constants of
% integration.

r_phi = sin(phi)*R;
lambda = sqrt(kappa/Sigma);

% first set up the matrix equation to find the constants of integration

M = [log(d/2/lambda), 1, besseli(0,d/2/lambda), besselk(0,d/2/lambda); ...
     2/d, 0, +besseli(1,d/2/lambda)/lambda, -besselk(1,d/2/lambda)/lambda; ...
    -1/r_phi^2, 0, 1/lambda^2*(besseli(0,r_phi/lambda)-1/r_phi*besseli(1,r_phi/lambda)), ...
            1/lambda^2*(besselk(0,r_phi/lambda)+1/r_phi*besselk(1,r_phi/lambda)); ...
     1/r_phi, 0, besseli(1,r_phi/lambda)/lambda, -besselk(1,r_phi/lambda)/lambda];

c = [0,0,1/R-sqrt(2*w/kappa),tan(phi)]';

C = M\c;

h = C(1)*log(r/lambda) + C(2) + C(3)*besseli(0,r/lambda) + C(4)*besselk(0,r/lambda);

A = 2*pi*trapz(r, ...
    r.*sqrt(1+(C(1)./r+C(3)/lambda*besseli(1,r/lambda)-C(4)/lambda*besselk(1,r/lambda)).^2))...
    + d^2*(1-pi/4);

% old equations

% M = [log(d/2), 1, besselj(0,d/2/lambda), bessely(0,d/2/lambda);...
%      2/d, 0, -besselj(1,d/2/lambda)/lambda, -bessely(1,d/2/lambda)/lambda;...
%      1/L0, 0, -besselj(1,L0/lambda)/lambda, -bessely(1,L0/lambda)/lambda;...
%      -1, 0, L0/lambda*(1+1/lambda^2)*besselj(1,L0/lambda), ...
%             L0/lambda*(1+1/lambda^2)*bessely(1,L0/lambda)];
% 
% c = [0,0,tan(phi),0]';
% 
% C = M\c;
% 
% h = C(1)*log(r) + C(2) + C(3)*besselj(0,r/lambda) + C(4)*bessely(0,r/lambda);
% 
% A = 2*pi*trapz(r, ...
%     r.*sqrt(1+(C(1)./r-C(3)/lambda*besselj(1,r/lambda)-C(4)/lambda*bessely(1,r/lambda)).^2))...
%     + d^2*(1-pi/4);

end