clear variables
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
R = 0.05; %um
d = 0.4; %um
%initial A/lipid = 0.58 nm^2 = 0.58e-6 um^2
A_lip_base = (0.58e-6);
N_lipid = d^2/A_lip_base;
% z = eps(1);
% Sigma_init = 0.7004;
% sigma = 0.399999975000000;
% sigma = 5;
% w = 4+sigma*2+0.05;
w = 1.2125*4;
% lambda_init = R/sqrt(Sigma_init);
kD = 300/10^12*1e9;     % picoJ/um^2
kappa = 1e-19*1e12;     % picoJ
epsilon = w*kappa/(2*R^2);
ws = epsilon-2*kappa/R^2;

% phi_vals = deg2rad(linspace(1, 179,50));
% phi_vals = 0.698131700797732;
phi_vals = deg2rad(170);
z_vals = 1-cos(phi_vals);
% z_vals = linspace(0+eps(1),1.9999, 50);
for jj = 1:length(z_vals)
z = z_vals(jj);
phi = acos(1-z);
% lambda = lambda_init;
% Sigma = Sigma_init;

A_init = pi*d^2 + 2*pi*R^2*z - pi*R^2*sin(phi)^2;
Sigma = get_Sigma(A_init, N_lipid, A_lip_base, kD);

while true
    % first find the shape of the membrane subject to some initial sigma;
    out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, psi_dot_init);
    
    solution = deval(out, linspace(0,out.xe, 1000));
    r_nonlin = solution(2,:);
    h_nonlin = solution(3,:)-solution(3,end);
    
    if out.x(end)<d/2
        r_nonlin = [solution(2,1:end-1),d/2];
        h_nonlin = [solution(3,1:end-1),solution(3,end)]-solution(3,end);
    elseif out.x(end)>d/2
        r_nonlin(r_nonlin>d/2) = d/2;
        h_nonlin(r_nonlin>d/2) = 0;
    end
    
    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = out.ye(8)*2*pi + pi*((d/2)^2-out.ye(2)^2) + d^2*(1-pi/4);
    
    A = S_A+S_B;

    % get the 
    
    E_all_nonlinear_self(1,ii) = E;
    E_all_nonlinear_self(2,ii) = E_adhesion;
    E_all_nonlinear_self(3,ii) = E_stretch_A;
    E_all_nonlinear_self(4,ii) = E_stretch_B;
    E_all_nonlinear_self(5,ii) = E_bend_A;
    E_all_nonlinear_self(6,ii) = E_bend_B;
    
    alpha_A_vals_nonlinear_self(ii) = alpha_A;
    alpha_B_vals_nonlinear_self(ii) = alpha_B;
    %     h_phi_vals_nonlinear_self(ii) = h_phi;
    S_A_vals_nonlinear_self(ii) = S_A;
    S_B_vals_nonlinear_self(ii) = S_B;
    Sigma_vals_nonlinear_self(ii) = Sigma;

end


end

function Sigma = get_Sigma(A, N_lipid, A_lip_base,kD)
    A_per_lipid = A/N_lipid;
    % assume tension is just linearly related to lipid density
    Sigma = (A_per_lipid-A_lip_base)*kD;
end