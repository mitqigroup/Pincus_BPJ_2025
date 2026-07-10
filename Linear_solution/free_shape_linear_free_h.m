function [C, delA, Ebend] = free_shape_linear_free_h(R, phi, kappa, Sigma)
% solves the shape equation for the free surface of a section of membrane
% bound at one end to a sphere of radius R at the point phi radians from
% the bottom of the sphere, with a flat surface at d/2. The shape is
% characterised by the length lambda, which is equal to the square root of
% the membrane bending energy kappa divided by the tension Sigma. Also
% outputs the total area of the membrane, as well as the constants of
% integration.

if Sigma > 0

    lambda = sqrt(kappa/Sigma);
    

    r_phi = R*sin(phi);
    
    C(1) = lambda*tan(phi)*besselk(0,r_phi/lambda,1)/besselk(1,r_phi/lambda,1);
    C(2) = -lambda*tan(phi)/besselk(1,r_phi/lambda);

    delA = -pi/2*r_phi*tan(phi)^2*(r_phi...
    -(besselk(0,r_phi/lambda,1)*(2*lambda*besselk(1,r_phi/lambda,1)...
    +r_phi*besselk(0,r_phi/lambda,1))/(besselk(1,r_phi/lambda,1)^2)));

    Ebend = 1/2*pi*r_phi^2*tan(phi)^2*kappa/lambda^2*(...
        -besselk(0,r_phi/lambda,1)^2+besselk(1,r_phi/lambda,1)^2)/besselk(1,r_phi/lambda,1)^2;

else 

    C(1) = 0;
    C(2) = 0;
    delA = 0;
    Ebend = 0;

end

