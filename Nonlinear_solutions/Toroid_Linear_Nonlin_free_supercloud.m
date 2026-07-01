clear variables
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');

% addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts");
% addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Linear_solution");

MyTaskID = 0;
NumberOfTasks = 1;

% check that the environment variables have been read in correctly
if ~(exist('MyTaskID', 'var')&&exist('NumberOfTasks', 'var'))
    error('Environment variables not set correctly')
end

% taskIDs count from zero, alter this here
MyTaskID = MyTaskID + 1;

% constants
R = 0.05;
sigma = 1e-4;
d = sqrt(pi*R^2/sigma);    % um
% d = 2;    % um
% sigma = pi*R^2/d^2;
phi = pi/12;
kD  = 300/10^12*1e9;    % picoJ/um^2
kappa = 1e-19*1e12;     % picoJ
kappa_bar = kappa/(kD*R^2);
alpha_i = 0.000316227766017;
% alpha_i = 0.0046;
% zeta = 0.00002;            % dimensionless
% epsilon = -zeta*kD;     % picoJ/um^2
% epsilon = -1e-2;
% epsilon = -0.000062*1.5;
% epsilon = -0.0784e-03;
% epsilon = -0.0809e-03;
epsilon = -1.501310728910000e-04;
% -0.0809
n0 = 1;                 % fraction
Sigma_bar = alpha_i*kD*R^2/kappa;
ep_bar = -epsilon*n0*R^2/(2*kappa);
% zeta = 0.98*2*kappa_bar         % dimensionless
% zeta = 1.001*(4*kappa_bar*(1+alpha_i))/(2-2*alpha_i-alpha_i^2+4*kappa_bar*(1+alpha_i))
% zeta = 5*(4*kappa_bar*(1+alpha_i))/(2-2*alpha_i-alpha_i^2+4*kappa_bar*(1+alpha_i))
% zeta = 2e-2
% zeta = 0.1*alpha_i^2/2
% gamma = 0.1*R;            % um
% epsilon = -kappa/gamma^2;
% kappa = 1e-17*1e12;     % picoJ
zeta = -epsilon*n0/kD;   % dimensionless
% alpha_i = sqrt(-zeta/2);        % fraction
% alpha_i = 0;
ten_param=-d*epsilon*n0/sqrt(kD*kappa);
psi_dot_init = -1;
N = 3e3;                % number of points in quadrature

% phi_vals = flip(deg2rad(linspace(0.1,60,48*4)));
phi_vals = deg2rad(linspace(5,175,4));

%% split up job between task IDs
size_phi_vals = size(phi_vals);
total_phi_vals = size_phi_vals(2);
min_sets_per_job = floor(total_phi_vals/NumberOfTasks);
extra_sets = mod(total_phi_vals, NumberOfTasks);
if MyTaskID <= extra_sets
    my_set_min = (min_sets_per_job+1)*(MyTaskID-1)+1
    my_set_max = (min_sets_per_job+1)*(MyTaskID)
else
    my_set_min = (min_sets_per_job+1)*(extra_sets) ...
        + min_sets_per_job*(MyTaskID-extra_sets-1) + 1
    my_set_max = (min_sets_per_job+1)*(extra_sets) ...
        + min_sets_per_job*(MyTaskID-extra_sets)
end
phi_vals_self = phi_vals(my_set_min:my_set_max);
my_set_size = size(phi_vals_self);
my_set_length = my_set_size(2);

save_name = sprintf('data/task%i_results.mat', MyTaskID-1);

%%
E_all = zeros(6,length(phi_vals_self));

E_all_toroid_self = nan(6,length(phi_vals_self));
alpha_A_vals_toroid_self = nan(1,length(phi_vals_self));
alpha_B_vals_toroid_self = nan(size(alpha_A_vals_toroid_self));
S_A_vals_toroid_self = nan(size(alpha_A_vals_toroid_self));
S_B_vals_toroid_self = nan(size(alpha_A_vals_toroid_self));
Sigma_vals_toroid_self = nan(size(alpha_A_vals_toroid_self));
rho_vals_self = nan(size(alpha_A_vals_toroid_self));

E_all_linear_self = nan(6,length(phi_vals_self));
alpha_A_vals_linear_self = nan(1,length(phi_vals_self));
alpha_B_vals_linear_self = nan(size(alpha_A_vals_toroid_self));
S_A_vals_linear_self = nan(size(alpha_A_vals_toroid_self));
S_B_vals_linear_self = nan(size(alpha_A_vals_toroid_self));
Sigma_vals_linear_self = nan(size(alpha_A_vals_toroid_self));

E_all_nonlinear_self = nan(6,length(phi_vals_self));
alpha_A_vals_nonlinear_self = nan(1,length(phi_vals_self));
alpha_B_vals_nonlinear_self = nan(size(alpha_A_vals_toroid_self));
h_phi_vals_nonlinear_self = nan(size(alpha_A_vals_toroid_self));
S_A_vals_nonlinear_self = nan(size(alpha_A_vals_toroid_self));
S_B_vals_nonlinear_self = nan(size(alpha_A_vals_toroid_self));
Sigma_vals_nonlinear_self = nan(size(alpha_A_vals_toroid_self));

tic
for ii = 1:length(phi_vals_self)
    phi = phi_vals_self(ii);
    rad2deg(phi)
    
    if ii==1
        alpha_B_init = alpha_i+0.01;
        alpha_A_init = (1+alpha_B_init)/(1+zeta)-1;
    else
        alpha_A_init = alpha_A;
        alpha_B_init = alpha_B;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linear solution, only for phi<85 degrees
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if phi<1.2
        % %% use initial stretch to solve for minimum attachment
        const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi];
            [out,~,~,~,~,~,~] = ...
                get_linear_minimum(const, [alpha_A_init, alpha_B_init]);
        
        alpha_A = out(1);
        alpha_B = out(2);
        
        % get the shape of the free region
        Sigma = kD*alpha_B;
        lambda = sqrt(kappa/Sigma);
        r_phi = sin(phi)*R;
        r = linspace(r_phi, d/2,N);
        
        [C, delA, Ebend] = free_shape_linear_free_h(R, phi, kappa, Sigma);
        
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = delA+(d^2-pi*r_phi^2);
        
        E_adhesion = epsilon*n0*S_A./(1+alpha_A);
        E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
        E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
        E_bend_B = Ebend;
        E_bend_A = 4*pi*kappa*(1-cos(phi));
        
        h = C(1)+C(2)*besselk(0,r/lambda);
        
        E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;
        
        E_all_linear_self(1,ii) = E;
        E_all_linear_self(2,ii) = E_adhesion;
        E_all_linear_self(3,ii) = E_stretch_A;
        E_all_linear_self(4,ii) = E_stretch_B;
        E_all_linear_self(5,ii) = E_bend_A;
        E_all_linear_self(6,ii) = E_bend_B;
        
        alpha_A_vals_linear_self(ii) = alpha_A;
        alpha_B_vals_linear_self(ii) = alpha_B;
        S_A_vals_linear_self(ii) = S_A;
        S_B_vals_linear_self(ii) = S_B;
        Sigma_vals_linear_self(ii) = Sigma;

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Toroidal solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e5, 'algorithm', 'sqp');   

    const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi];
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        get_toroid_minimum(const, [alpha_A_init, alpha_B_init, R/20]);
    
    alpha_A = out(1);
    alpha_B = out(2);
    rho = out(3);
    
    delta = sin(phi)*(R+rho);
    
    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = d^2 - pi*delta^2 + 2*pi*rho*(delta*phi-rho*(1-cos(phi)));
    
    E_adhesion = epsilon*S_A./(1+alpha_A);
    E_stretch_A = kD/2*alpha_A.^2*S_A./(1+alpha_A);
    E_stretch_B = kD/2*alpha_B.^2*S_B./(1+alpha_B);
    E_bend_A = 4*pi*kappa*(1-cos(phi));
    a = delta/rho;
    if a>1
        E_bend_B = 2*pi*kappa*(a^2/sqrt(a^2-1)*...
            (acot(sqrt(a^2-1))+atan((a*tan(phi/2)-1)/sqrt(a^2-1)))-2*(1-cos(phi)));
    else
        E_bend_B = 2*pi*kappa*(a^2/sqrt(1-a^2)*...
            (-acoth(sqrt(1-a^2))+atanh((1-a*tan(phi/2))/sqrt(1-a^2)))-2*(1-cos(phi)));
    end
    
    E = E_adhesion+E_stretch_A+E_stretch_B+E_bend_A+E_bend_B;
    
    E_all_toroid_self(1,ii) = E;
    E_all_toroid_self(2,ii) = E_adhesion;
    E_all_toroid_self(3,ii) = E_stretch_A;
    E_all_toroid_self(4,ii) = E_stretch_B;
    E_all_toroid_self(5,ii) = E_bend_A;
    E_all_toroid_self(6,ii) = E_bend_B;
    
%     h_phi = sin(-3*pi/2)*rho+rho*cos(phi)- (sin(-3*pi/2+phi)*rho+rho*cos(phi));
    
    alpha_A_vals_toroid_self(ii) = alpha_A;
    alpha_B_vals_toroid_self(ii) = alpha_B;
%     h_phi_vals_toroid(ii) = h_phi;
    S_A_vals_toroid_self(ii) = S_A;
    S_B_vals_toroid_self(ii) = S_B;
    Sigma_vals_toroid_self(ii) = alpha_B*kD;
    rho_vals_self(ii) = rho;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nonlinear solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi];
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        get_nonlinear_minimum(const, [alpha_A_init, alpha_B_init]);
    
    alpha_A = out(1);
    alpha_B = out(2);
    
    lambda_stretch = lam_vals.eqnonlin;
    
    % get the shape of the free region
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
    
    out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, false);

    solution = deval(out, linspace(0,out.x(end), 1000));
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
    S_B = out.y(8,end)*2*pi + pi*((d/2)^2-out.y(2,end)^2) + d^2*(1-pi/4);

    E_adhesion = epsilon*n0*S_A./(1+alpha_A);
    E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
    E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
    E_bend_B = out.y(7,end)*kappa*pi;
    E_bend_A = 4*pi*kappa*(1-cos(phi));
    E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;
    
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
toc

MyTaskID = MyTaskID-1;
save(save_name);

%% functions

function [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_linear_minimum(constants, inputs)

    epsilon = constants(1);
    n0 = constants(2);
    d = constants(3);
    R = constants(4);
    kD = constants(5);
    kappa = constants(6);
    alpha_i = constants(7);
    N = constants(8);
    phi = constants(9);

    grad = [];
    hessian = [];
    
    alpha_A_init = inputs(1);
    alpha_B_init = inputs(2);

    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;
    
    Sigma = kD*alpha_B;
    r_phi = sin(phi)*R;
    lambda = sqrt(kappa/Sigma);
    lamt = lambda*sqrt(alpha_B);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = (d^2-pi*r_phi^2);
    S_i = d^2;

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%         'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
%         'SpecifyObjectiveGradient', true, 'CheckGradients', true);

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%         'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
%         'SpecifyObjectiveGradient', true);

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%         'OptimalityTolerance', 1e-14, 'ConstraintTolerance', 1e-14, 'StepTolerance', 1e-14,...
%         'SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient', true,...
%         'Diagnostics', 'on', 'Display', 'iter-detailed');

    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 4e3, 'algorithm', 'sqp',...
        'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
        'SpecifyObjectiveGradient', false);

    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init],...
        [],[],[],[],[-0.1,-0.000001],[0.1, 0.1], ...
        @constraint, options);
% 
%     clear options
%     options.verify_level = 1;
%     options.printfile = 'snsolve_print.txt';
%     options.scale_option = 2;
%     options.major_feasibility_tolerance = 1e-12;
%     options.major_optimality_tolerance = 1e-12;
% %     options.function_precision = 1e-8;
% %     options.major_optimality_tolerance = 1e-4;
%     [out,fval,exitflag,output,lam_vals,states] = ...
%         snsolve(@objective,[alpha_A_init, alpha_B_init], [],[],[],[],...
%         [-0.1,-0.1]',[0.1, 0.1]', ...
%         @constraint, options);

    function [f,g] = objective(y_obj)

        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        
        Sigma = kD*alpha_B;
        r_phi = sin(phi)*R;
        lambda = sqrt(kappa/Sigma);
        lamt = lambda*sqrt(alpha_B);

        [~, delA, Ebend] =...
            free_shape_linear_free_h(R, phi, kappa, Sigma);
    
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = delA+(d^2-pi*r_phi^2);

        alpha_A_gradient = (kD/2*S_A*alpha_A*(alpha_A+2)-epsilon*n0*S_A)...
            /(alpha_A+1)^2;

%         alpha_A_gradient = -epsilon*n0*S_A./(1+alpha_A)^2+kD*alpha_A*S_A/(1+alpha_A)...
%             -kD*1/2*alpha_A^2*S_A/(1+alpha_A)^2

        parts(1) = -kD*1/2*d^2*alpha_B^2/(1+alpha_B)^2;
        parts(2) = kD*d^2*alpha_B/(1+alpha_B);
        parts(3) = kD/2*pi*r_phi^2*alpha_B^2/(1+alpha_B)^2;
        parts(4) = -kD*pi*r_phi^2*alpha_B/(1+alpha_B);
        parts(5) = -kD*pi*r_phi^2*alpha_B^2*tan(phi)^2/(1+alpha_B)^2/4;
        parts(6) = +kD*3/4*r_phi^2*alpha_B*tan(phi)^2/(1+alpha_B)^2;
        parts(7) = -pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2);
        parts(8) = +(pi*besselk(0,r_phi/lambda,1)* ...
        (-4*kD*r_phi*sqrt(alpha_B)*lamt^4 ...
         -2*kD*r_phi*alpha_B^(3/2)*lamt^4 ...
         -2*r_phi^3*sqrt(alpha_B)*kappa ...
         -4*r_phi^3*alpha_B^(3/2)*kappa ...
         -2*r_phi^3*alpha_B^(5/2)*kappa ...
         +kD*r_phi^3 *alpha_B^(3/2)*lamt^2 ...
         +kD*r_phi^3*alpha_B^(5/2)*lamt^2)*tan(phi)^2)...
         /(4*(1+alpha_B)^2*lamt^3*besselk(1,r_phi/lambda,1));
        parts(9) = (besselk(0,r_phi/lambda,1)^2)* ...
         ((kD*pi*r_phi^2*tan(phi)^2*alpha_B^2)/(4*(1+alpha_B)^2) ...
         -(kD*5*pi*r_phi^2*tan(phi)^2*alpha_B)/(8*(1+alpha_B)) ...
         +pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2) ...
         +(pi*r_phi^2*(4*kappa+4*kD*alpha_B*kappa-3*kD*alpha_B*lamt^2)*tan(phi)^2 ...
         /(8*(1+alpha_B)*lamt^2))...
         )/(besselk(1,r_phi/lambda,1))^2;
        parts(10) = (besselk(0,r_phi/lambda,1)^3* ...
         (pi*r_phi^3*sqrt(alpha_B)*kappa*tan(phi)^2/(4*lamt^3)  ...
         -kD*pi*r_phi^3*alpha_B^(3/2)*tan(phi)^2/(8*(1+alpha_B)*lamt) ...
         +pi*r_phi^3*tan(phi)^2*(2*sqrt(alpha_B)*kappa+2*alpha_B^(3/2) ...
         *kappa-kD*alpha_B^(3/2)*lamt^2)/(8*(1+alpha_B)*lamt^3))) ...
         /(besselk(1,r_phi/lambda,1)^3);

        alpha_B_gradient = (...
        -kD*1/2*d^2*alpha_B^2/(1+alpha_B)^2 ...
        +kD*d^2*alpha_B/(1+alpha_B) ...
        +kD/2*pi*r_phi^2*alpha_B^2/(1+alpha_B)^2 ... 
        -kD*pi*r_phi^2*alpha_B/(1+alpha_B) ...
        -kD*pi*r_phi^2*alpha_B^2*tan(phi)^2/(1+alpha_B)^2/4 ...
        +kD*3/4*r_phi^2*alpha_B*tan(phi)^2/(1+alpha_B)^2 ...
        -pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2)...
        +(pi*besselk(0,r_phi/lambda,1)* ...
        (-4*kD*r_phi*sqrt(alpha_B)*lamt^4 ...
         -2*kD*r_phi*alpha_B^(3/2)*lamt^4 ...
         -2*r_phi^3*sqrt(alpha_B)*kappa ...
         -4*r_phi^3*alpha_B^(3/2)*kappa ...
         -2*r_phi^3*alpha_B^(5/2)*kappa ...
         +kD*r_phi^3 *alpha_B^(3/2)*lamt^2 ...
         +kD*r_phi^3*alpha_B^(5/2)*lamt^2)*tan(phi)^2)...
         /(4*(1+alpha_B)^2*lamt^3*besselk(1,r_phi/lambda,1)) ...
        + (besselk(0,r_phi/lambda,1)^2)* ...
         ((kD*pi*r_phi^2*tan(phi)^2*alpha_B^2)/(4*(1+alpha_B)^2) ...
         -(kD*5*pi*r_phi^2*tan(phi)^2*alpha_B)/(8*(1+alpha_B)) ...
         +pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2) ...
         +(pi*r_phi^2*(4*kappa+4*kD*alpha_B*kappa-3*kD*alpha_B*lamt^2)*tan(phi)^2 ...
         /(8*(1+alpha_B)*lamt^2))...
         )/(besselk(1,r_phi/lambda,1))^2 ...
        + (besselk(0,r_phi/lambda,1)^3* ...
         (pi*r_phi^3*sqrt(alpha_B)*kappa*tan(phi)^2/(4*lamt^3)  ...
         -kD*pi*r_phi^3*alpha_B^(3/2)*tan(phi)^2/(8*(1+alpha_B)*lamt) ...
         +pi*r_phi^3*tan(phi)^2*(2*sqrt(alpha_B)*kappa+2*alpha_B^(3/2) ...
         *kappa-kD*alpha_B^(3/2)*lamt^2)/(8*(1+alpha_B)*lamt^3))) ...
         /(besselk(1,r_phi/lambda,1)^3));

        g = [alpha_A_gradient, alpha_B_gradient]';

        E_adhesion = epsilon*n0*S_A./(1+alpha_A);
        E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
        E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
        E_bend_B = Ebend;
        E_bend_A = 4*pi*kappa*(1-cos(phi));

        % stretching, adhesion and bending energy
        f = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;

    end

    function [c,ceq, gc, gceq]= constraint(~)
    
        c = [];
        ceq = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

        gc = [];

%         alpha_A
%         alpha_B
%         S_A
%         S_B
%         ceq

%         parts(1) = 

        alpha_A_gradient = -S_A/(1+alpha_A)^2;
%         alpha_B_gradient_old = 1/(1+alpha_B)^2*( ...
%             -d^2+pi*r_phi^2-pi*r_phi^2*tan(phi)^2/2 ...
%             +(pi*besselk(0,r_phi/lambda,1))/(2*sqrt(alpha_B)*lamt*besselk(1,r_phi/lambda,1)) ...
%             *(2*r_phi*lamt^2+r_phi^3+r_phi^3*alpha_B)*tan(phi)^2) ...
%             + 1/(1+alpha_B)*(pi*r_phi^2*tan(phi)^2/2 ...
%             -pi*r_phi^3*tan(phi)^2*besselk(0,r_phi/lambda,1)^3 ...
%             /besselk(1,r_phi/lambda,1)^3/(2*sqrt(alpha_B))) ...
%             +besselk(0,r_phi/lambda,1)^2/besselk(1,r_phi/lambda,1)^2 ...
%             *(pi*r_phi^2*tan(phi)/(2*(1+alpha_B)^2)-pi*r_phi^2*tan(phi)/(alpha_B*(1+alpha_B)));

        alpha_B_gradient = (-1).*d.^2.*(1+alpha_B).^(-2)+pi.*R.^2.*(1+alpha_B).^(-2).*sin(phi).^2+(-1/2) ...
          .*pi.*R.^2.*(1+alpha_B).^(-2).*sin(phi).^2.*tan(phi).^2+(1/2).*pi.*R.^2.* ...
          alpha_B.^(-1).*(1+alpha_B).^(-1).*sin(phi).^2.*tan(phi).^2+(-1/2).*pi.*R.^3.* ...
          alpha_B.^(-1/2).*(1+alpha_B).^(-1).*lamt.^(-1).*besselk(0,R.*alpha_B.^(1/2).*lamt.^( ...
          -1).*sin(phi),1).^3.*besselk(1,R.*alpha_B.^(1/2).*lamt.^(-1).*sin(phi),1).^(-3).* ...
          sin(phi).^3.*tan(phi).^2+(1/2).*pi.*alpha_B.^(-1/2).*(1+alpha_B).^(-2).*lamt.^(-1) ...
          .*besselk(0,R.*alpha_B.^(1/2).*lamt.^(-1).*sin(phi),1).*besselk(1,R.*alpha_B.^( ...
          1/2).*lamt.^(-1).*sin(phi),1).^(-1).*(2.*R.*lamt.^2.*sin(phi)+R.^3.*sin(phi) ...
          .^3+R.^3.*alpha_B.*sin(phi).^3).*tan(phi).^2+besselk(0,R.*alpha_B.^(1/2).*lamt.^( ...
          -1).*sin(phi),1).^2.*besselk(1,R.*alpha_B.^(1/2).*lamt.^(-1).*sin(phi),1).^(-2).* ...
          ((1/2).*pi.*R.^2.*(1+alpha_B).^(-2).*sin(phi).^2.*tan(phi).^2+(-1).*pi.* ...
          R.^2.*alpha_B.^(-1).*(1+alpha_B).^(-1).*sin(phi).^2.*tan(phi).^2);


        gceq = [alpha_A_gradient,  alpha_B_gradient]';


    end

end

function [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_nonlinear_minimum(constants, inputs)

    epsilon = constants(1);
    n0 = constants(2);
    d = constants(3);
    R = constants(4);
    kD = constants(5);
    kappa = constants(6);
    alpha_i = constants(7);
    N = constants(8);
    phi = constants(9);
    
    alpha_A_init = inputs(1);
    alpha_B_init = inputs(2);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = 0;
    S_i = d^2;

    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;

    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
        'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12);
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init],...
        [],[],[],[],[-0.1,eps(1)],[0.1, 0.1], ...
        @constraint, options);

    function f = objective(y_obj)

        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        
        Sigma = kD*alpha_B;
        lambda = sqrt(kappa/Sigma);

        psi_dot_init = -sqrt(2*(1-cos(phi)))/lambda;
        
        out_shape = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, false);
        
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = out_shape.y(8,end)*2*pi + pi*((d/2)^2-out_shape.y(2,end)^2) + d^2*(1-pi/4);
        
        % stretching, adhesion and bending energy
        f = epsilon*n0*S_A./(1+alpha_A) ...
          + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
          + out_shape.y(7,end)*kappa*pi + 4*pi*kappa*(1-cos(phi));

    end

    function [c,ceq]= constraint(~)
    
        c(1) = 0;
        ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
%         c(2) = 0;
%         ceq(2) = alpha_B-alpha_A;
    
    end

end

function [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_toroid_minimum(constants, inputs)

    epsilon = constants(1);
    n0 = constants(2);
    d = constants(3);
    R = constants(4);
    kD = constants(5);
    kappa = constants(6);
    alpha_i = constants(7);
    N = constants(8);
    phi = constants(9);
    
    alpha_A_init = inputs(1);
    alpha_B_init = inputs(2);
    rho_init = inputs(3);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = 0;
    S_i = d^2;

    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;
    rho = rho_init;

    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
        'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12);
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init, rho_init],...
        [],[],[],[],[-0.1,-0.1, 0],[0.1, 0.1, 5*d], ...
        @constraint, options);

    function f = objective(y_obj)
        
        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        rho = y_obj(3); %toroid radius
    
        delta = sin(phi)*(R+rho);
    
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = d^2 - pi*delta^2 + 2*pi*rho*(delta*phi-rho*(1-cos(phi)));
    
        E_adhesion = epsilon*S_A./(1+alpha_A);
        E_stretch_A = kD/2*alpha_A.^2*S_A./(1+alpha_A);
        E_stretch_B = kD/2*alpha_B.^2*S_B./(1+alpha_B);
        E_bend_A = 4*pi*kappa*(1-cos(phi));
        a = delta/rho;
        if a>1
            E_bend_B = 2*pi*kappa*(a^2/sqrt(a^2-1)*...
                (acot(sqrt(a^2-1))+atan((a*tan(phi/2)-1)/sqrt(a^2-1)))-2*(1-cos(phi)));
        else
            E_bend_B = 2*pi*kappa*(a^2/sqrt(1-a^2)*...
                (-acoth(sqrt(1-a^2))+atanh((1-a*tan(phi/2))/sqrt(1-a^2)))-2*(1-cos(phi)));
        end
    
        f = E_adhesion+E_stretch_A+E_stretch_B+E_bend_A+E_bend_B;

    end

    function [c,ceq]= constraint(~)
    
        c(1) = 0;
        ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
    
    end

end 
