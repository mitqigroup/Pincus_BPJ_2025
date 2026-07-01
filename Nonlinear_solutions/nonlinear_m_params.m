clear variables
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');

% addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts");
% addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Linear_solution");
% addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Nonlinear_solutions");

MyTaskID = 0;
NumberOfTasks = 1;

% check that the environment variables have been read in correctly
if ~(exist('MyTaskID', 'var')&&exist('NumberOfTasks', 'var'))
    error('Environment variables not set correctly')
end

% taskIDs count from zero, alter this here
MyTaskID = MyTaskID + 1;

% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = [0.05];                 % um
sigma_vals = 0.00021;                % surface fraction
%d_vals = sqrt(R_vals.^2./sigma_vals);           % um
kD_vals = 0.3;                        % picoJ/um^2
zeta_vals = 0.001;                  % dimensionless
epsilon_vals = -zeta_vals*kD_vals(1);              % picoJ/um^2
n0_vals = 1;                                    % fraction
kappa_vals = 1e-7;         % picoJ
% kappa_vals = logspace(-21,-15, 60)*1e12;
alpha_i_vals = [1e-5];                            % fraction
m_vals = [0];
% phi_vals = flip(deg2rad(linspace(0.01,179.9,150)));
phi_vals = deg2rad(linspace(1,179,6));
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
                            for mm=1:length(m_vals)
                            jj = jj+1;
                                parameter_set(jj, 1:8) = [epsilon_vals(ee),...
                                    n0_vals(nn),...
                                    sqrt(pi*R_vals(rr)^2/sigma_vals(ss)),...
                                    R_vals(rr),...
                                    kD_vals(kk),...
                                    kappa_vals(pp),...
                                    alpha_i_vals(aa),...
                                    m_vals(mm)];
                            end
                        end
                    end
                end
            end
        end
    end
end

% randomly shuffle the parameter set so that we have a more even
% distribution between each core on our node. We can unshuffle it later on
% if we want to. We MUST set the random number generator to the same value,
% since it must be identical between cores
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

E_all_nonlinear_self = nan(6,length(phi_vals_self), my_set_length);
alpha_A_vals_nonlinear_self = nan(length(phi_vals_self), my_set_length);
alpha_B_vals_nonlinear_self = nan(size(alpha_A_vals_nonlinear_self));
h_phi_vals_nonlinear_self = nan(size(alpha_A_vals_nonlinear_self));
S_A_vals_nonlinear_self = nan(size(alpha_A_vals_nonlinear_self));
S_B_vals_nonlinear_self = nan(size(alpha_A_vals_nonlinear_self));
Sigma_vals_nonlinear_self = nan(size(alpha_A_vals_nonlinear_self));
phi_end_nonlinear_self = nan(size(alpha_A_vals_nonlinear_self));

for ii = 1:my_set_length
    epsilon = parameter_set_self(ii, 1);
    n0 = parameter_set_self(ii, 2);
    d = parameter_set_self(ii, 3);
    R = parameter_set_self(ii, 4);
    kD = parameter_set_self(ii, 5);
    kappa = parameter_set_self(ii, 6);
    alpha_i = parameter_set_self(ii, 7);
    m = parameter_set_self(ii, 8);
    fprintf('parameter set is %i \n', ii);
    fprintf('epsilon = %0.4g \n', epsilon);
    fprintf('n0 = %0.4g \n', n0);
    fprintf('d = %0.4g \n', d);
    fprintf('R = %0.4g \n', R);
    fprintf('kD = %0.4g \n', kD);
    fprintf('kappa = %0.4g \n', kappa);
    fprintf('alpha_i = %0.4g \n', alpha_i);
    fprintf('m = %0.4g \n', m);
    zeta = epsilon*n0/kD;
    sigma = sqrt(pi*R^2/d^2);

    tic
    for jj = 1:length(phi_vals_self)
        phi = phi_vals_self(jj)
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
    
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Nonlinear solution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi, m];
        [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
            get_nonlinear_minimum(const, [alpha_A_init, alpha_B_init]);
        
        alpha_A = out(1);
        alpha_B = out(2);
        
        lambda_stretch = lam_vals.eqnonlin;
        
        % get the shape of the free region
        Sigma = kD*alpha_B;
        lambda = sqrt(kappa/Sigma);
        r_phi = sin(phi)*R;
        
        out = Canalejo_2021_m_func(R, phi, Sigma, kappa, m, d, false);
    
        solution = deval(out, linspace(0,out.x(end), 1000));
        r_nonlin = solution(2,:);
        h_nonlin = solution(4,:)-solution(4,end);
    
        if out.x(end)<d/2
            r_nonlin = [solution(2,1:end-1),d/2];
            h_nonlin = [solution(3,1:end-1),solution(3,end)]-solution(3,end);
        elseif out.x(end)>d/2
            r_nonlin(r_nonlin>d/2) = d/2;
            h_nonlin(r_nonlin>d/2) = 0;
        end
        
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = out.y(6,end) + d^2 - pi*R^2*sin(phi)^2;
        r_end = out.ye(2);
    
        E_adhesion = epsilon*n0*S_A./(1+alpha_A);
        E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
        E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
        E_bend_B = out.y(7,end) - 2*pi*kappa*(r_end^2-R^2*sin(phi)^2)*m^2;
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
    m = constants(10);
    
    alpha_A_init = inputs(1);
    alpha_B_init = inputs(2);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = 0;
    S_i = d^2;

    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;

    options = optimoptions('fmincon','MaxFunEvals', 1e3, 'MaxIter',...
        1e3, 'algorithm', 'sqp',...
        'OptimalityTolerance', 1e-12, 'ConstraintTolerance',...
        1e-12, 'StepTolerance', 1e-12, 'Display', 'iter');
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init],...
        [],[],[],[],[-10,eps(1)],[10, 10], ...
        @constraint, options);

    function f = objective(y_obj)

        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        
        Sigma = kD*alpha_B;
        
%         out_shape = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, false);
        out_shape = Canalejo_2021_m_func(R, phi, Sigma, kappa, m, d, false);

        r_end = out_shape.ye(2);
        
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = out_shape.y(6,end) + d^2 - pi*R^2*sin(phi)^2;
        
        % stretching, adhesion and bending energy
        f = epsilon*n0*S_A./(1+alpha_A) ...
          + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
          + out_shape.y(7,end) - 2*pi*kappa*(r_end^2-R^2*sin(phi)^2)*m^2 ...
          + 4*pi*kappa*(1-cos(phi));

    end

    function [c,ceq]= constraint(~)
    
        c(1) = 0;
        ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
%         c(2) = 0;
%         ceq(2) = alpha_B-alpha_A;
    
    end

end
