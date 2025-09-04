function [coeff, error, Q] = pfit2(p, x, y, dy)
    %Implements polyfit scheme in 15.4 Numerical Methods 3rd Ed
        % p is horizontal array containing polynomial powers
        % e.g. [4 3 2 1 0] for 4th order polynomial
        % x is timestep widths (or general x-data)
        % y is dependent data
        % dy is standard error in each y point
    
    V = bsxfun(@power, x, p);
    A = V./dy;
    b = y./dy;
    alpha = A'*A;
    beta = A'*b;
    coeff = alpha^-1*beta;
    error = sqrt(diag(alpha^-1));
    
    pi = V*coeff;
    Q = sum((y-pi).^2./dy.^2);
end