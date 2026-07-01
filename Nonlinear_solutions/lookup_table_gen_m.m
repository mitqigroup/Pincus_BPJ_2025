% clear variables
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

% actually everything depends only on Sigma and phi, but we'll write out
% some example values just to be safe
R = 1;
kappa = 1e-7;
d = 1000;
sigma_bar_vals = logspace(-5,8,16);
phi_vals = smoothstep(linspace(0.03,1-0.03,20),3)*pi;
m_bar_vals = [-logspace(1,-3, 23), 0, logspace(-3,1, 24)];

%% split up job between task IDs
size_m = size(m_bar_vals);
total_m = size_m(2);
min_sets_per_job = floor(total_m/NumberOfTasks);
extra_sets = mod(total_m, NumberOfTasks);
if MyTaskID <= extra_sets
    my_set_min = (min_sets_per_job+1)*(MyTaskID-1)+1;
    my_set_max = (min_sets_per_job+1)*(MyTaskID);
else
    my_set_min = (min_sets_per_job+1)*(extra_sets) ...
        + min_sets_per_job*(MyTaskID-extra_sets-1) + 1;
    my_set_max = (min_sets_per_job+1)*(extra_sets) ...
        + min_sets_per_job*(MyTaskID-extra_sets);
end
m_self = m_bar_vals(my_set_min:my_set_max);
my_set_size = size(m_self);
my_set_length = my_set_size(1);

save_name = sprintf('data/task%i_results.mat', MyTaskID-1);

%% run
lm = length(m_self);
lp = length(phi_vals);
ls = length(sigma_bar_vals);

delA_bar = nan(lm,lp,ls);
E_bend_bar = nan(size(delA_bar));
rend_bar = nan(size(delA_bar));
psi_end = nan(size(delA_bar));

tic
% ticBytes(gcp);
% parfor (dd = 1:ld, 8)
for mm = 1:lm
    m_bar = m_self(mm)*R;
    for ss = 1:ls
        sigma_bar = sigma_bar_vals(ss);
        if mod(ss,1)==0
            fprintf('%i \n', ss);
        end
        for pp = 1:lp
            phi = phi_vals(pp);
            fprintf('%0.4g \n', phi);
            Sigma = sigma_bar*kappa/R^2;
            m = m_bar/R;
            out = Canalejo_2021_m_func(R, phi, Sigma, kappa, m, d, false);
            r = out.y(2,:);
            rend = r(end);
            rphi = r(1);    
            delA_bar(mm,pp,ss) = out.y(6,end)/R^2;
            E_bend_bar(mm,pp,ss) = out.y(7,end)/kappa ...
                - 2*pi*(rend^2-rphi^2)*m^2;
            rend_bar(mm,pp,ss) = rend/R;
            psi_end(mm,pp,ss) = out.y(1,end);
        end
    end
end
% tocBytes(gcp)
toc

results_self.delA_bar = delA_bar;
results_self.E_bend_bar = E_bend_bar;
results_self.rend_bar = rend_bar;
results_self.psi_end = psi_end;

results_self_mat(1,:,:,:) = delA_bar;
results_self_mat(2,:,:,:) = E_bend_bar;
results_self_mat(3,:,:,:) = rend_bar;
results_self_mat(4,:,:,:) = psi_end;

save(save_name, 'results_self_mat', "results_self");

