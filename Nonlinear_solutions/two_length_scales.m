% clear variables
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');

% lookup_table = load("lookup_table.mat")
% lookup_table = load("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Nonlinear_solutions/gridded_interp.mat")
% lookup_table = load("gridded_interp.mat");

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

data = load("ExpDATA.mat");

expData = data.ExpDATA;

% 
% R = 0.1;
% d = 10;
% kD = 0.3;
% epsilon = -1e-4;
% kappa = 1e-7;
% alpha_i = 1e-4;
% Sig_b = 1e-5;
% zeta = -epsilon/kD;
% phi_vals = deg2rad(linspace(1,179, 50));
% alpha_A = 0.02;
% alpha_B = 0.01;

kB = 1.380649e-23;
NA = 6.0221408e+23;
T = 298.15;
l = 16e-9;
M = 1e5;
kappa = 33*kB*T;
Sigma = 1e-9;
c = 0.04;
n = expData(:,2)/100/M*NA*1e6;
kD = 2e-5;

omega = n*l*kB*T;
omega_star = n*kB*T.*(l-3*((c*kB*T)./(4*n*kappa)).^(1/3));
Rp = expData(:,1);
Rv = mean(expData(:,3:4),2);

w_bar = 2*omega.*Rp.^2/(kappa);
w_bar_star = 2*omega_star.*Rp.^2/(kappa);
Sigma_bar = Sigma*Rp.^2/kappa;
type = expData(:,5);
sigma = 1/4*Rp.^2./Rv.^2;
Sigma_bar_end = (Sigma+sigma*4*0.3).*Rp.^2/kappa;

% assuming vesicles are prolate spheroids
c_s = expData(:,3);
a_s = expData(:,4);
e_s = sqrt(1-a_s.^2./c_s.^2);
A_sp = 2*pi*a_s.^2.*(1+c_s./(a_s.*e_s).*asin(e_s));
V_sp = 4/3*pi*a_s.^2.*c_s;
R_vsp = sqrt(A_sp/(4*pi));
R_vol_vsp = (3*V_sp/(4*pi)).^(1/3);
red_V_sp = 3*V_sp./(4*pi*R_vsp.^3);
excess_A_sp = (A_sp-4*pi*R_vol_vsp.^2)./A_sp;

% assuming vesicles are oblate spheroids
a_s = expData(:,3);
c_s = expData(:,4);
e_s = sqrt(1-c_s.^2./a_s.^2);
A_so = 2*pi*a_s.^2.*(1+(1-e_s.^2)./e_s.*atanh(e_s));
V_so = 4/3*pi*a_s.^2.*c_s;
R_vso = sqrt(A_so/(4*pi));
R_vol_vso = (3*V_so/(4*pi)).^(1/3);
red_V_so = 3*V_so./(4*pi*R_vso.^3);
excess_A_so = (A_so-4*pi*R_vol_vso.^2)./A_so;

parameter_set = [];
for Sig_b = logspace(-10,-6,20)
    parameter_set_1 = [-omega_star,...
        ones(size(omega_star)),...
        sqrt(A_so)*1e6,...
        Rp*1e6,...
        ones(size(omega_star))*kD,...
        ones(size(omega_star))*kappa*1e12,...
        ones(size(omega_star))*1e-10,...
        ones(size(omega_star))*Sig_b];
    
    parameter_set = [parameter_set; parameter_set_1];
end

phi_vals = deg2rad(linspace(0.1,179.9,200));
%phi_vals = [phi_vals(70:-1:1),phi_vals(71:length(phi_vals))];
N = 1e3;

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

for ii=1:length(parameter_set_self)
    epsilon = parameter_set_self(ii, 1);
    n0 = parameter_set_self(ii, 2);
    d = parameter_set_self(ii, 3);
    R = parameter_set_self(ii, 4);
    kD = parameter_set_self(ii, 5);
    kappa = parameter_set_self(ii, 6);
    alpha_i = parameter_set_self(ii, 7);
    Sig_b = parameter_set_self(ii, 8);
    zeta = -epsilon/kD;
    alpha_A = 0.02;
    alpha_B = 0.01;
    tic
    for pp = 1:length(phi_vals)
        phi = phi_vals(pp);
        % first solve the deserno problem
        out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sig_b, false);
        r = out.y(2,:);
        h = out.y(3,:);
        rend = r(end);
        rphi = r(1);
        delA = 2*pi*out.y(8,end) - pi*(rend.^2-rphi.^2);
        E_bend_free = out.y(7,end)*kappa*pi;
    %     E(pp) = E_bend_free;
    %     delA = 0;
    %     E_bend_free = 0;
        
        % now solve the microplastics problem
        options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
            'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
            'SpecifyObjectiveGradient', false,...
            'SpecifyConstraintGradient', false,'Diagnostics', 'off',...
            'Display', 'none');
    
        const = [zeta, alpha_i, d, R, phi, delA];
        [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
            fmincon(@(y) free_phi(y, const),[alpha_A, alpha_B], [],[],[],[],...
            [-0.1,-0.1],[0.1, 0.1], ...
            @(y) area_con_phi(y,const), options);
        
        alpha_A = out(1);
        alpha_B = out(2);
        
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = d^2 - pi*R^2*sin(phi)^2;
        
        E_self(ii,pp) = fval*kD+E_bend_free+4*pi*kappa*(1-cos(phi));
        aA_self(ii,pp) = alpha_A;
        aB_self(ii,pp) = alpha_B;
    end
    toc
%     save(save_name);
end

% figure();
% plot(phi_vals/pi, E);

%% functions
function [f, g] = free_phi(y, const)
    % function of free energy to minimise

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);
    phi = const(5);
    delA = const(6);

    alpha_A = y(1);
    alpha_B = y(2);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = d^2 - pi*R^2*sin(phi)^2 + delA;

    g = [-zeta*S_A./(1+alpha_A)^2+alpha_A*S_A/(1+alpha_A)-1/2*alpha_A^2*S_A/(1+alpha_A)^2;
         alpha_B*S_B/(1+alpha_B)-1/2*alpha_B^2*S_B/(1+alpha_B)^2];

%         parts(1) = -kD*1/2*d^2*alpha_B^2/(1+alpha_B)^2;
%         parts(2) = kD*d^2*alpha_B/(1+alpha_B);
%         parts(3) = kD/2*pi*r_phi^2*alpha_B^2/(1+alpha_B)^2;
%         parts(4) = -kD*pi*r_phi^2*alpha_B/(1+alpha_B);

    f = -zeta*S_A./(1+alpha_A) ...
      + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
      - alpha_i^2*d^2/(1+alpha_i)/2;

end

function [c,ceq, gc, gceq] = area_con_phi(y, const)

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);
    phi = const(5);
    delA = const(6);

    alpha_A = y(1);
    alpha_B = y(2);

    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = d^2 - pi*R^2*sin(phi)^2 + delA;

    c = [];

    ceq = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

    gc = [];

    gceq = [-S_A/(1+alpha_A)^2  -S_B/(1+alpha_B)^2]';

end