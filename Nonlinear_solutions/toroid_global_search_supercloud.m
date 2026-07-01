clear variables
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');

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
% sigma_vals = logspace(-4,-1,6);                % surface fraction
sigma_vals = logspace(-4,-1,6);                % surface fraction
%d_vals = sqrt(pi*R_vals.^2./sigma_vals);           % um
kD_vals = 300/10^12*1e9;                        % picoJ/um^2
kD_base = 300/10^12*1e9;                        % picoJ/um^2
% zeta_vals = logspace(-6,-5,1);                  % dimensionless
zeta_vals = 1e-3;                  % dimensionless
epsilon_vals = -zeta_vals*kD_base;              % picoJ/um^2
n0_vals = 1;                                    % fraction
kappa_vals = 1e-19*1e12;         % picoJ
% kappa_vals = logspace(-21,-15, 60)*1e12;
alpha_i_vals = 0.003;                            % fraction
% other constants
N = 3e5;

ii = 0;
for rr = 1:length(R_vals)
    for ss = 1:length(sigma_vals)
        for kk = 1:length(kD_vals)
            for ee = 1:length(epsilon_vals)
                for nn = 1:length(n0_vals)
                    for pp = 1:length(kappa_vals)
                        for aa = 1:length(alpha_i_vals)
                            ii = ii+1;
                            parameter_set(ii, 1:7) = [epsilon_vals(ee),...
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

%% run batch job for this task set
for ii = 1:my_set_length
    epsilon = parameter_set_self(ii, 1);
    n0 = parameter_set_self(ii, 2);
    d = parameter_set_self(ii, 3);
    R = parameter_set_self(ii, 4);
    kD = parameter_set_self(ii, 5);
    kappa = parameter_set_self(ii, 6);
    alpha_i = parameter_set_self(ii, 7);
    zeta = epsilon*n0/kD;
    sigma = sqrt(R^2/d^2);

    % solve microplastics for initial stretch of A and B
    options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e4);
    const = [zeta, alpha_i, d, R];
    [out,fval,~,~,lam_vals,grad,hessian] = ...
        fmincon(@(y) free_phi(y, const),[0.2, 0.2, pi/3], [],[],[],[],...
        [-1,-1,0],[Inf, Inf, pi/2], ...
        @(y) lipid_con_phi(y,const), options);
    alpha_A_init = out(1);
    alpha_B_init = out(2);
    phi_init = out(3);
    S_A_init = 2*pi*R^2*(1-cos(phi_init));
    S_B_init = d^2 - pi*R^2*sin(phi_init)^2;

    E_adhesion_init = epsilon*n0*S_A_init ./(1+alpha_A_init );
    E_stretch_A_init  = kD/2*(alpha_A_init .^2*S_A_init ./(1+alpha_A_init ));
    E_stretch_B_init  = kD/2*(alpha_B_init .^2*S_B_init ./(1+alpha_B_init ));
    E_total_init = E_adhesion_init+E_stretch_A_init+E_stretch_B_init;

    constants = [epsilon, n0, d, R, kD, kappa, alpha_i];
    inputs = [alpha_A_init, alpha_B_init, phi_init, R];
   
    [out,fval2,exitflag,output,solutions] = ...
        get_toroid_minimum(constants, inputs);

    solution_set{ii} = solutions;
    
    % get energies etc
    alpha_A = out(1);
    alpha_B = out(2);
    phi = out(3);
    rho = out(4);

    delta = sin(phi)*(R+rho);
    Sigma = kD*alpha_B;
    r_phi = sin(phi)*R;
    lambda = sqrt(kappa/Sigma);
    
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
    
    E_all(1,ii) = E;
    E_all(2,ii) = E_adhesion;
    E_all(3,ii) = E_stretch_A;
    E_all(4,ii) = E_stretch_B;
    E_all(5,ii) = E_bend_A;
    E_all(6,ii) = E_bend_B;
    
    alpha_A_vals(ii) = alpha_A;
    alpha_B_vals(ii) = alpha_B;
    phi_vals(ii) = phi;
    rho_vals(ii) = rho;
    S_A_vals(ii) = S_A;
    S_B_vals(ii) = S_B;
    Sigma_vals(ii) = Sigma;

    % do one round of the full nonlinear calculation with this solution,
    % and store the results for later checking.

    out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, 0);

    S_B_check(ii) = out.y(8)*2*pi + pi*((d/2)^2-out.y(2)^2) + d^2*(1-pi/4);
    E_bend_check(ii) = out.y(7,end)*pi*kappa;
    
end

save(save_name)

%% extra functions
function f = free_phi(y, const)
    % function of free energy to minimise

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);

    alpha_A = y(1);
    alpha_B = y(2);
    phi = y(3);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = d^2 - pi*R^2*sin(phi)^2;

    f = zeta*S_A./(1+alpha_A) ...
      + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B));

end

function [c,ceq] = lipid_con_phi(y, const)

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);

    alpha_A = y(1);
    alpha_B = y(2);
    theta = y(3);

    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(theta));
    S_B = d^2 - pi*R^2*sin(theta)^2;

    c(1) = 0;

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end

function [out,fval2,exitflag,output,solutions] = ...
    get_toroid_minimum(constants, inputs)

    epsilon = constants(1);
    n0 = constants(2);
    d = constants(3);
    R = constants(4);
    kD = constants(5);
    kappa = constants(6);
    alpha_i = constants(7);
    
    alpha_A_init = inputs(1);
    alpha_B_init = inputs(2);
    phi_init = inputs(3);
    rho_init = inputs(4);

    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;
    phi = phi_init;
    rho = rho_init;

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = 0;
    S_i = d^2;

    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
        'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
        'SpecifyObjectiveGradient', false);

    problem = createOptimProblem("fmincon",...
        'x0', [alpha_A_init, alpha_B_init, phi_init, rho_init],...
        'objective', @objective,...
        'lb', [-0.1,-0.1,0, 0],...
        'ub', [0.1, 0.1, pi, 5*d],...
        'nonlcon', @constraint,...
        'options', options);
    alpha_A_vals_inp = [-0.1, -alpha_i, -1e-3, 0, 1e-3, alpha_i, 0.1];
    phi_vals_inp = deg2rad([0, 5, 20, 45, 60, 85, 90,120,145,170,180]);
    rho_vals_inp = linspace(0,10*d,10);
%     alpha_A_vals_inp = [0, 1e-3, alpha_i, 0.1];
%     phi_vals_inp = deg2rad([0,45,90]);
    ptmatrix_this_run =...
            zeros(length(alpha_A_vals_inp),length(phi_vals_inp), 3);
    tic
    for aa = 1:length(alpha_A_vals_inp)
        for pp = 1:length(phi_vals_inp)
            for rr = 1:length(rho_vals_inp)
                phi = phi_vals_inp(pp);
                rho = rho_vals_inp(rr);
                alpha_B_vals_inp(pp,aa,rr) =...
                    fzero(@lipid_con_bend_test_A, rand()*0.2-0.1);
                ptmatrix_this_run(aa,pp,rr,1:4) = ...
                    [alpha_A_vals_inp(aa), alpha_B_vals_inp(pp,aa,rr),...
                     phi_vals_inp(pp), rho_vals_inp(rr)];
            end
        end
    end
    toc
    ptmatrix = reshape(ptmatrix_this_run, [numel(ptmatrix_this_run)/4, 4]);
    tpoints = CustomStartPointSet(ptmatrix);
    rs = RandomStartPointSet('NumStartPoints',100);
    gs = MultiStart("FunctionTolerance",1e-3, "XTolerance", 1e-3);
    tic
    [out,fval2,exitflag,output,solutions] = run(gs,problem, {tpoints,rs});
    toc

    function f = objective(y_obj)
        
        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        phi = y_obj(3);
        rho = y_obj(4); %toroid radius
    
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

    function ceq = lipid_con_bend_test_A(alpha_B)
   
        % problems when phi = 0, so change slightly
        if phi==0
            phi = phi+eps(1);
        end
        
        delta = sin(phi)*(R+rho);
    
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = d^2 - pi*delta^2 + 2*pi*rho*(delta*phi-rho*(1-cos(phi)));
    
        ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
    
    end

end
