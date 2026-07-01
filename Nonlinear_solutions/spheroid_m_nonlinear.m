%clear variables
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');

% lookup_table = load("lookup_table.mat")
lookup_table = load("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Nonlinear_solutions/gridded_interp.mat");
% lookup_table = load("gridded_interp.mat");

addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts");
addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Linear_solution");
addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Nonlinear_solutions");

%MyTaskID = 0;
%NumberOfTasks = 1;

% check that the environment variables have been read in correctly
if ~(exist('MyTaskID', 'var')&&exist('NumberOfTasks', 'var'))
    error('Environment variables not set correctly')
end

% taskIDs count from zero, alter this here
MyTaskID = MyTaskID + 1;

% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = [0.1];                 % um
sigma_vals = [1e-10,1e-5,5e-2];                % surface fraction
%d_vals = sqrt(R_vals.^2./sigma_vals);           % um
% d = 3.1623;
kD_vals = logspace(-2,2,8);                        % picoJ/um^2
% zeta_vals = 1e-5;                  % dimensionless
epsilon_vals = -logspace(-9,-2,32);              % picoJ/um^2
n0_vals = 1;                                    % fraction
kappa_vals = logspace(-9,-5,8);         % picoJ
% kappa_vals = logspace(-21,-15, 60)*1e12;
alpha_i_vals = [1e-7];                            % fraction
% phi_vals = flip(deg2rad(linspace(0.01,179.9,5)));
phi_vals = deg2rad(linspace(0.1,179.9,200));
%phi_vals = [phi_vals(70:-1:1),phi_vals(71:length(phi_vals))];
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
            alpha_B_init = alpha_i+sigma;
            alpha_A_init = (1+alpha_B_init)/(1+zeta)-1;
        else
            alpha_A_init = alpha_A;
            alpha_B_init = alpha_B;
        end

%         assert(alpha_A_init<0.1);
%         assert(alpha_B_init<0.1);
%         assert(alpha_B_init>0.99*alpha_i);

        fprintf('alpha_A init = %0.4g  \t', alpha_A_init);
        fprintf('alpha_B init = %0.4g \n', alpha_B_init);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Nonlinear solution lookup table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if phi>pi/2&&jj~=1
            % check if this is going to lead to an unphysical result. If it is,
            % let's just skip and move on
            Sigma = alpha_B_init*kD;
            out = free_shape_m_ellipse(a,b, d, phi, kappa, Sigma, m, false);
            psi_end = out.y(1,end);
            if psi_end>0.1
                % now we just skip
                E_all_nonlinear_self(:,jj,ii) = E_all_nonlinear_self(:,jj-1,ii);
                alpha_A_vals_nonlinear_self(jj,ii) = alpha_A_vals_nonlinear_self(jj-1,ii);
                alpha_B_vals_nonlinear_self(jj,ii) = alpha_B_vals_nonlinear_self(jj-1,ii);
                S_A_vals_nonlinear_self(jj,ii) = S_A_vals_nonlinear_self(jj-1,ii);
                S_B_vals_nonlinear_self(jj,ii) = S_B_vals_nonlinear_self(jj-1,ii);
                Sigma_vals_nonlinear_self(jj,ii) = Sigma_vals_nonlinear_self(jj-1,ii);
                phi_end_nonlinear_self(jj,ii) = psi_end;
                continue
            end
        end

        const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi];
        [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
            get_nonlinear_minimum(const, [alpha_A_init, alpha_B_init], lookup_table);

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
    
    end
    save(save_name, '-regexp', '^(?!(lookup_table)$).');
    toc
end

MyTaskID = MyTaskID-1;
save(save_name, '-regexp', '^(?!(lookup_table)$).');

%% functions

function [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_nonlinear_minimum(constants, inputs)

    epsilon = constants(1);
    n0 = constants(2);
    d = constants(3);
    a = constants(4);
    kD = constants(5);
    kappa = constants(6);
    alpha_i = constants(7);
    N = constants(8);
    phi = constants(9);
    c = constants(10);
    m = constants(11);

    beta_vals = atan(a/c*tan(linspace(0,phi,1000)));
    beta_vals(beta_vals<0) = beta_vals(beta_vals<0)+pi;
    t1 = (c^2-a^2)*cos(beta_vals).^2;
    integrand = (c*(2*a^2+t1)./(2*a.*(a^2+t1).^(3/2))*2-m).^2.*...
        a.*sin(beta_vals).*sqrt(a^2*cos(beta_vals).^2+c^2.*sin(beta_vals).^2);
    % t1 = (a^2-c^2)*cos(2*beta_vals);
    % integrand = (c*(3*a^2+c^2+t1)./(sqrt(2)*a.*(a^2+c^2+t1).^(3/2))-m).^2.*...
    %     a.*sin(beta_vals).*sqrt(a^2*cos(beta_vals).^2+c^2.*sin(beta_vals).^2);
    E_A = pi*kappa*trapz(beta_vals, integrand);
    S_A = 2*pi*trapz(beta_vals, a.*sin(beta_vals).*sqrt(a^2*cos(beta_vals).^2+c^2.*sin(beta_vals).^2));
    
    alpha_A_init = inputs(1);
    alpha_B_init = inputs(2);

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
        
        out_shape = free_shape_m_ellipse(a, c, d, phi, kappa, Sigma, m, false);
        
        S_B = out_shape.y(8,end)*2*pi + pi*((d/2)^2-out_shape.y(2,end)^2) + d^2*(1-pi/4);
        
        % stretching, adhesion and bending energy
        f = epsilon*n0*S_A./(1+alpha_A) ...
          + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
          + out_shape.y(7,end)*kappa*pi + E_A;

    end

    function [c,ceq]= constraint(~)
    
        c(1) = 0;
        ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
%         c(2) = 0;
%         ceq(2) = alpha_B-alpha_A;
    
    end

end