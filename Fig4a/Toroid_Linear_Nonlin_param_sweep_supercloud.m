% clear variables
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');

addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts");
addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Linear_solution");
addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Nonlinear_solutions");

% MyTaskID = 0;
% NumberOfTasks = 1;

% check that the environment variables have been read in correctly
if ~(exist('MyTaskID', 'var')&&exist('NumberOfTasks', 'var'))
    error('Environment variables not set correctly')
end

% taskIDs count from zero, alter this here
MyTaskID = MyTaskID + 1;

% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = [0.05];                 % um
sigma_vals = [logspace(-9,log10(0.00013),10),linspace(0.00018,0.00026,28),logspace(log10(0.00036),-1,10)];
%d_vals = sqrt(R_vals.^2./sigma_vals);           % um
kD_vals = 300/10^12*1e9;                        % picoJ/um^2
kD_base = 300/10^12*1e9;                        % picoJ/um^2
zeta_vals = 0.001;                  % dimensionless
epsilon_vals = -zeta_vals*kD_base;              % picoJ/um^2
n0_vals = 1;                                    % fraction
kappa_vals = 1e-19*1e12;         % picoJ
% kappa_vals = logspace(-21,-15, 60)*1e12;
alpha_i_vals = [1e-5];                            % fraction
phi_vals = flip(deg2rad(linspace(0.01,179.9,150)));
% other constants
N = 3e5;

jj = 0;
for rr = 1:length(R_vals)
    for ss = 1:length(sigma_vals)
        for kk = 1:length(kD_vals)
            for ee = 1:length(epsilon_vals)
                for nn = 1:length(n0_vals)
                    for pp = 1:length(kappa_vals)
                        for aa = 1:length(alpha_i_vals)
                            jj = jj+1;
                            parameter_set(jj, 1:7) = [epsilon_vals(ee),...
                                n0_vals(nn),...
                                sqrt(pi*R_vals(rr)^2/sigma_vals(ss)),...
                                R_vals(rr),...
                                kD_vals(kk),...
                                kappa_vals(pp),...
                                alpha_i_vals(aa)];
                        end
                    end
                end
            end
        end
    end
end

% randomly shuffle the parameter set so that we have a more even distribution between
% each core on our node. We can unshuffle it later on if we want to. We MUST set the
% random number generator to the same value, since it must be identical between cores
rng('default');
rng(1);

size_params = size(parameter_set);

permutation_array = randperm(size_params(1));

parameter_set_permuted = parameter_set;
% parameter_set_unpermuted = parameter_set;

for ii=1:size_params(2)
    parameter_set_permuted(:,ii) = parameter_set(permutation_array,ii);
end
parameter_set_original = parameter_set;
parameter_set = parameter_set_permuted;

% for ii=1:size_params(2)
    % for jj = 1:length(permutation_array)
        % kk = permutation_array(jj);
        % parameter_set_unpermuted(kk,ii) = parameter_set_permuted(jj,ii);
    % end
% end

%% split up job between task IDs
size_parameter_set = size(parameter_set);
total_parameter_sets = size_parameter_set(1);
min_sets_per_job = floor(total_parameter_sets/NumberOfTasks);
extra_sets = mod(total_parameter_sets, NumberOfTasks);
if MyTaskID <= extra_sets
    my_set_min = (min_sets_per_job+1)*(MyTaskID-1)+1
    my_set_max = (min_sets_per_job+1)*(MyTaskID)
else
    my_set_min = (min_sets_per_job+1)*(extra_sets) ...
        + min_sets_per_job*(MyTaskID-extra_sets-1) + 1
    my_set_max = (min_sets_per_job+1)*(extra_sets) ...
        + min_sets_per_job*(MyTaskID-extra_sets)
end
parameter_set_self = parameter_set(my_set_min:my_set_max, :);
my_set_size = size(parameter_set_self);
my_set_length = my_set_size(1);

save_name = sprintf('data/task%i_results.mat', MyTaskID-1);

%%

phi_vals_self = phi_vals;

E_all = zeros(6,length(phi_vals_self), my_set_length);

E_all_toroid_self = nan(6,length(phi_vals_self), my_set_length);
alpha_A_vals_toroid_self = nan(length(phi_vals_self), my_set_length);
alpha_B_vals_toroid_self = nan(size(alpha_A_vals_toroid_self));
S_A_vals_toroid_self = nan(size(alpha_A_vals_toroid_self));
S_B_vals_toroid_self = nan(size(alpha_A_vals_toroid_self));
Sigma_vals_toroid_self = nan(size(alpha_A_vals_toroid_self));
rho_vals_self = nan(size(alpha_A_vals_toroid_self));

E_all_linear_self = nan(6,length(phi_vals_self), my_set_length);
alpha_A_vals_linear_self = nan(length(phi_vals_self), my_set_length);
alpha_B_vals_linear_self = nan(size(alpha_A_vals_toroid_self));
S_A_vals_linear_self = nan(size(alpha_A_vals_toroid_self));
S_B_vals_linear_self = nan(size(alpha_A_vals_toroid_self));
Sigma_vals_linear_self = nan(size(alpha_A_vals_toroid_self));

E_all_nonlinear_self = nan(6,length(phi_vals_self), my_set_length);
alpha_A_vals_nonlinear_self = nan(length(phi_vals_self), my_set_length);
alpha_B_vals_nonlinear_self = nan(size(alpha_A_vals_toroid_self));
h_phi_vals_nonlinear_self = nan(size(alpha_A_vals_toroid_self));
S_A_vals_nonlinear_self = nan(size(alpha_A_vals_toroid_self));
S_B_vals_nonlinear_self = nan(size(alpha_A_vals_toroid_self));
Sigma_vals_nonlinear_self = nan(size(alpha_A_vals_toroid_self));
phi_end_nonlinear_self = nan(size(alpha_A_vals_toroid_self));

for ii = 1:my_set_length
    epsilon = parameter_set_self(ii, 1);
    n0 = parameter_set_self(ii, 2);
    d = parameter_set_self(ii, 3);
    R = parameter_set_self(ii, 4);
    kD = parameter_set_self(ii, 5);
    kappa = parameter_set_self(ii, 6);
    alpha_i = parameter_set_self(ii, 7);
    fprintf('parameter set is %i \n', ii);
    fprintf('epsilon = %0.4g \n', epsilon);
    fprintf('n0 = %0.4g \n', n0);
    fprintf('d = %0.4g \n', d);
    fprintf('R = %0.4g \n', R);
    fprintf('kD = %0.4g \n', kD);
    fprintf('kappa = %0.4g \n', kappa);
    fprintf('alpha_i = %0.4g \n', alpha_i);
    zeta = epsilon*n0/kD;
    sigma = sqrt(pi*R^2/d^2);

    tic
    for jj = 1:length(phi_vals_self)
        phi = phi_vals_self(jj);
        fprintf('phi degreees = %0.4g \n', rad2deg(phi));
        
        if jj==1
            alpha_B_init = alpha_i+0.01;
            alpha_A_init = (1+alpha_B_init)/(1+zeta)-1;
        else
            alpha_A_init = alpha_A;
            alpha_B_init = alpha_B;
        end

        fprintf('alpha_A init = %0.4g  \t', alpha_A_init);
        fprintf('alpha_B init = %0.4g \n', alpha_B_init);
        
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
            
            E_all_linear_self(1,jj,ii) = E;
            E_all_linear_self(2,jj,ii) = E_adhesion;
            E_all_linear_self(3,jj,ii) = E_stretch_A;
            E_all_linear_self(4,jj,ii) = E_stretch_B;
            E_all_linear_self(5,jj,ii) = E_bend_A;
            E_all_linear_self(6,jj,ii) = E_bend_B;
            
            alpha_A_vals_linear_self(jj,ii) = alpha_A;
            alpha_B_vals_linear_self(jj,ii) = alpha_B;
            S_A_vals_linear_self(jj,ii) = S_A;
            S_B_vals_linear_self(jj,ii) = S_B;
            Sigma_vals_linear_self(jj,ii) = Sigma;
    
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Toroidal solution, doesn't work for very very high phi
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if phi<3.09
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
            
            E_all_toroid_self(1,jj,ii) = E;
            E_all_toroid_self(2,jj,ii) = E_adhesion;
            E_all_toroid_self(3,jj,ii) = E_stretch_A;
            E_all_toroid_self(4,jj,ii) = E_stretch_B;
            E_all_toroid_self(5,jj,ii) = E_bend_A;
            E_all_toroid_self(6,jj,ii) = E_bend_B;
            
        %     h_phi = sin(-3*pi/2)*rho+rho*cos(phi)- (sin(-3*pi/2+phi)*rho+rho*cos(phi));
            
            alpha_A_vals_toroid_self(jj,ii) = alpha_A;
            alpha_B_vals_toroid_self(jj,ii) = alpha_B;
        %     h_phi_vals_toroid(jj,ii) = h_phi;
            S_A_vals_toroid_self(jj,ii) = S_A;
            S_B_vals_toroid_self(jj,ii) = S_B;
            Sigma_vals_toroid_self(jj,ii) = alpha_B*kD;
            rho_vals_self(jj,ii) = rho;

        end
    
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
        
        E_all_nonlinear_self(1,jj,ii) = E;
        E_all_nonlinear_self(2,jj,ii) = E_adhesion;
        E_all_nonlinear_self(3,jj,ii) = E_stretch_A;
        E_all_nonlinear_self(4,jj,ii) = E_stretch_B;
        E_all_nonlinear_self(5,jj,ii) = E_bend_A;
        E_all_nonlinear_self(6,jj,ii) = E_bend_B;
        
        alpha_A_vals_nonlinear_self(jj,ii) = alpha_A;
        alpha_B_vals_nonlinear_self(jj,ii) = alpha_B;
    %     h_phi_vals_nonlinear_self(jj,ii) = h_phi;
        S_A_vals_nonlinear_self(jj,ii) = S_A;
        S_B_vals_nonlinear_self(jj,ii) = S_B;
        Sigma_vals_nonlinear_self(jj,ii) = Sigma;
        phi_end_nonlinear_self(jj,ii) = out.y(1,end);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
    save(save_name);
    toc
end

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
        'SpecifyObjectiveGradient', false, 'Display', 'off');

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
        'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12, 'Display', 'off');
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init],...
        [],[],[],[],[-0.1,eps(1)],[0.1, 0.1], ...
        @constraint, options);

    function f = objective(y_obj)

        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        
        Sigma = kD*alpha_B;
        
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
        'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12, 'Display', 'off');
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
