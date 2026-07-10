function [C, A, E, E_bend] = free_shape_linear_no_curve(r_phi, d, phi, kappa, Sigma, h_phi)
% solves the shape equation for the free surface of a section of membrane
% bound at one end to a sphere of radius R at the point phi radians from
% the bottom of the sphere, with a flat surface at d/2. The shape is
% characterised by the length lambda, which is equal to the square root of
% the membrane bending energy kappa divided by the tension Sigma. Also
% outputs the total area of the membrane, as well as the constants of
% integration.

if Sigma==0
    Sigma = Sigma+1e-10;
end

if Sigma > 0

    lambda = sqrt(kappa/Sigma);

    if d/2/lambda>100

        % asymptotic solution
        C4_bar = (d/(2*r_phi*lambda)*besseli(1,d/2/lambda,1)*h_phi ...
            - besseli(0,d/2/lambda,1)*tan(phi) -...
            log(2*r_phi/d)*d/(2*lambda)*besseli(1,d/2/lambda,1)*tan(phi))/...
            (d/(2*r_phi*lambda)*besselk(0,r_phi/lambda,1)*besseli(1,d/2/lambda,1)...
            + 1/lambda*besselk(1,r_phi/lambda,1)*besseli(0,d/2/lambda,1)...
            + log(2*r_phi/d)*d/(2*lambda^2)*besselk(1,r_phi/lambda,1)*besseli(1,d/2/lambda,1));
        C3_bar = (tan(phi)+1/lambda*besselk(1,r_phi/lambda,1)*C4_bar)/...
            (-d/(2*r_phi*lambda)*besseli(1,d/2/lambda,1));

        area_func = @(x) x.*sqrt(1+(-besseli(1,d/2/lambda,1)*C3_bar*d/(2*lambda)*1./x...
            +exp((x-d/2)/lambda).*C3_bar/lambda.*besseli(1,x/lambda,1)...
            -exp((r_phi-x)/lambda).*C4_bar/lambda.*besselk(1,x/lambda,1)).^2);
        bend_func = @(x) x.*(exp((x-d/2)/lambda).*C3_bar/lambda^2.*besseli(0,x/lambda,1)...
            +exp((r_phi-x)/lambda).*C4_bar/lambda^2.*besselk(0,x/lambda,1)).^2;
        sig_func = @(x) x.*(-besseli(1,d/2/lambda,1)*C3_bar*d/(2*lambda)*1./x...
            +exp((x-d/2)/lambda).*C3_bar/lambda.*besseli(1,x/lambda,1)...
            -exp((r_phi-x)/lambda).*C4_bar/lambda.*besselk(1,x/lambda,1)).^2;
        
        % these are wrong! We actually don't want to solve them directly in
        % the kappa -> 0 limit, since we need asymptotics
        C(1) = 0;
        C(2) = 0;
        C(3) = h_phi/log(2*r_phi/d);
        C(4) = -h_phi*log(2*r_phi/lambda)/log(d/2/lambda);

        AbsTol = 1e-7;
        RelTol = 1e-4;
        
        %solve for area and energy
        E_bend = kappa/2*2*pi*integral(bend_func,r_phi,d/2,...
            'AbsTol',AbsTol,'RelTol',RelTol);
        E = E_bend + Sigma/2*2*pi*integral(sig_func,r_phi,d/2,...
            'AbsTol',AbsTol,'RelTol',RelTol);
        A = 2*pi*integral(area_func,r_phi,d/2,...
            'AbsTol',AbsTol,'RelTol',RelTol) + d^2*(1-pi/4);

    else
        A_c = d/(2*r_phi*lambda)*besselk(1,d/2/lambda) - 1/lambda*besselk(1,r_phi/lambda);
        B_c = -d/(2*r_phi*lambda)*besseli(1,d/2/lambda) + 1/lambda*besseli(1,r_phi/lambda);
        D_c = besseli(0,r_phi/lambda) - besseli(0,d/2/lambda) ...
            - log(2*r_phi/d)*d/2/lambda*besseli(1,d/2/lambda);
        E_c = besselk(0,r_phi/lambda) - besselk(0,d/2/lambda) ...
            + log(2*r_phi/d)*d/2/lambda*besselk(1,d/2/lambda);
        
        C(4) = (h_phi - D_c*tan(phi)/B_c)/(-A_c*D_c/B_c+E_c);
        C(3) = (tan(phi) - C(4)*A_c)/B_c;
        C(1) = -C(3)*d/2/lambda*besseli(1,d/2/lambda) + C(4)*d/2/lambda*besselk(1,d/2/lambda);
        C(2) = -C(1)*log(d/2/lambda) - C(3)*besseli(0,d/2/lambda) - C(4)*besselk(0, d/2/lambda);

        area_func = @(x) x.*sqrt(1+(C(1)./x+C(3)/lambda*besseli(1,x/lambda)-C(4)/lambda*besselk(1,x/lambda)).^2);
        bend_func = @(x) x.*(C(3)/lambda^2*besseli(0,x/lambda)+C(4)/lambda^2*besselk(0,x/lambda)).^2;
        sig_func = @(x) x.*(C(1)./x+C(3)/lambda*besseli(1,x/lambda)-C(4)/lambda*besselk(1,x/lambda)).^2;

        AbsTol = 1e-7;
        RelTol = 1e-4;
        
        %solve for area and energy
        E_bend = kappa/2*2*pi*integral(bend_func,r_phi,d/2,...
            'AbsTol',AbsTol,'RelTol',RelTol);
        E = E_bend + Sigma/2*2*pi*integral(sig_func,r_phi,d/2,...
            'AbsTol',AbsTol,'RelTol',RelTol);
        A = 2*pi*integral(area_func,r_phi,d/2,...
            'AbsTol',AbsTol,'RelTol',RelTol) + d^2*(1-pi/4);

    end

elseif Sigma < 0

    lambda = sqrt(kappa/-Sigma);
    
    H_c = d/(2*r_phi*lambda)*besselj(1,d/2/lambda) - 1/lambda*besselj(1,r_phi/lambda);
    I_c = d/(2*r_phi*lambda)*bessely(1,d/2/lambda) - 1/lambda*bessely(1,r_phi/lambda);
    F_c = besselj(0,r_phi/lambda) - besselj(0,d/2/lambda) ...
        + log(2*r_phi/d)*d/2/lambda*besselj(1,d/2/lambda);
    G_c = bessely(0,r_phi/lambda) - bessely(0,d/2/lambda) ...
        + log(2*r_phi/d)*d/2/lambda*bessely(1,d/2/lambda);
    
    C(4) = (F_c*tan(phi)-H_c*h_phi)/(I_c*F_c-G_c*H_c);
    C(3) = (h_phi - C(4)*G_c)/F_c;
    C(1) = C(3)*d/2/lambda*besselj(1,d/2/lambda) + C(4)*d/2/lambda*bessely(1,d/2/lambda);
    C(2) = -C(1)*log(d/2/lambda) - C(3)*besselj(0,d/2/lambda) - C(4)*bessely(0, d/2/lambda);

    if lambda>5e-2

        area_func = @(x) x.*sqrt(1+(C(1)./x-C(3)/lambda*besselj(1,x/lambda)-C(4)/lambda*bessely(1,x/lambda)).^2);
        bend_func = @(x) x.*(-C(3)/lambda^2*besselj(0,x/lambda)-C(4)/lambda^2*bessely(0,x/lambda)).^2;
        sig_func = @(x) x.*(C(1)./x-C(3)/lambda*besselj(1,x/lambda)-C(4)/lambda*bessely(1,x/lambda)).^2;

        AbsTol = 1e-7;
        RelTol = 1e-4;
        
        %solve for area and energy
        E_bend = kappa/2*2*pi*integral(bend_func,r_phi,d/2,...
            'AbsTol',AbsTol,'RelTol',RelTol);
        E = E_bend + Sigma/2*2*pi*integral(sig_func,r_phi,d/2,...
            'AbsTol',AbsTol,'RelTol',RelTol);
        A = 2*pi*integral(area_func,r_phi,d/2,...
            'AbsTol',AbsTol,'RelTol',RelTol) + d^2*(1-pi/4);
    else
        N = 1e4;
        r = linspace(r_phi, d/2, N);

        hderiv = C(1)./r-C(3)/lambda*besselj(1,r/lambda)-C(4)/lambda*bessely(1,r/lambda);
        lap_h = -C(3)/lambda^2*besselj(0,r/lambda)-C(4)/lambda^2*bessely(0,r/lambda);
        
        A = 2*pi*trapz(r, r.*sqrt(1+(hderiv).^2)) + d^2*(1-pi/4);
        E_bend = kappa/2*2*pi*trapz(r, r.*lap_h.^2);
        E = E_bend + Sigma/2*2*pi*trapz(r, r.*hderiv.^2);
    end

end

% AbsTol = 1e-7;
% RelTol = 1e-4;
% 
% %solve for area and energy
% E_bend = kappa/2*2*pi*integral(bend_func,r_phi,d/2,...
%     'AbsTol',AbsTol,'RelTol',RelTol);
% E = E_bend + Sigma/2*2*pi*integral(sig_func,r_phi,d/2,...
%     'AbsTol',AbsTol,'RelTol',RelTol);
% A = 2*pi*integral(area_func,r_phi,d/2,...
%     'AbsTol',AbsTol,'RelTol',RelTol) + d^2*(1-pi/4);



