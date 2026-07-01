% clear variables
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');

% addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts");
% addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Linear_solution");
% addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts/Nonlinear_solutions");

MyTaskID = 4;
NumberOfTasks = 48;

% check that the environment variables have been read in correctly
if ~(exist('MyTaskID', 'var')&&exist('NumberOfTasks', 'var'))
    error('Environment variables not set correctly')
end

% taskIDs count from zero, alter this here
MyTaskID = MyTaskID + 1;

% actually everything depends only on Sigma and phi, but we'll write out
% some example values just to be safe
alpha_i = 1e-4;
R = 1;
kappa = 1e-7;
sigma_bar_vals = logspace(-5,8,800);
% phi_vals = linspace(0.01,pi-0.1,100);
%phi_vals = deg2rad(linspace(1,179,10));
phi_vals = smoothstep(linspace(0.03,1-0.03,200),3)*pi;
d_bar_vals = logspace(log10(3*R), 5,48);
d_bar_vals = d_bar_vals(2:end);
% d = 0.4;

%% split up job between task IDs
size_d = size(d_bar_vals);
total_d = size_d(2);
min_sets_per_job = floor(total_d/NumberOfTasks);
extra_sets = mod(total_d, NumberOfTasks);
if MyTaskID <= extra_sets
    my_set_min = (min_sets_per_job+1)*(MyTaskID-1)+1
    my_set_max = (min_sets_per_job+1)*(MyTaskID)
else
    my_set_min = (min_sets_per_job+1)*(extra_sets) ...
        + min_sets_per_job*(MyTaskID-extra_sets-1) + 1
    my_set_max = (min_sets_per_job+1)*(extra_sets) ...
        + min_sets_per_job*(MyTaskID-extra_sets)
end
d_self = d_bar_vals(my_set_min:my_set_max);
my_set_size = size(d_self);
my_set_length = my_set_size(1);

save_name = sprintf('data/task%i_results.mat', MyTaskID-1);


%% run
ld = length(d_self);
lp = length(phi_vals);
ls = length(sigma_bar_vals);

delA_bar = nan(ld,lp,ls);
E_bend_bar = nan(size(delA_bar));
rend_bar = nan(size(delA_bar));
phi_end = nan(size(delA_bar));

tic
% ticBytes(gcp);
% parfor (dd = 1:ld, 8)
for dd = 1:ld
    d = d_self(dd)*R;
    for ss = 600:602
        sigma_bar = sigma_bar_vals(ss);
        if mod(ss,1)==0
            fprintf('%i \n', ss);
        end
        for pp = 1:lp
            phi = phi_vals(pp);
            fprintf('%0.4g \n', phi);
            Sigma = sigma_bar*kappa/R^2;
            out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, false);
            r = out.y(2,:);
            rend = r(end);
            rphi = r(1);
            A = out.y(8,end);
            phi_end(dd,pp,ss) = out.y(1,end);        
            delA_bar(dd,pp,ss) = (2*pi*A - pi*(rend.^2-rphi.^2))/R^2;
            E_bend_bar(dd,pp,ss) = out.y(7,end);
            rend_bar(dd,pp,ss) = rend/R;
        end
    end
end
% tocBytes(gcp)
toc

results_self.phi_end = phi_end;
results_self.delA_bar = delA_bar;
results_self.E_bend_bar = E_bend_bar;
results_self.rend_bar = rend_bar;

results_self_mat(1,:,:,:) = phi_end;
results_self_mat(2,:,:,:) = delA_bar;
results_self_mat(3,:,:,:) = E_bend_bar;
results_self_mat(4,:,:,:) = rend_bar;

save(save_name, 'results_self_mat', "results_self");

%% combine variables back together
num_tasks = 48;
results = [];
for ii=1:num_tasks-1
    load_name = sprintf('data/task%i_results.mat', ii );
    data = load(load_name);
    results = cat(2,results,data.results_self_mat);
end

phi_end = squeeze(results(1,:,:,:));
delA_bar = squeeze(results(2,:,:,:));
E_bend_bar = squeeze(results(3,:,:,:));
rend_bar = squeeze(results(4,:,:,:));

%% plotting results
% % figure();
% hold on
% axes1 = gca;
% [sigma_mesh,phi_mesh] = meshgrid(sigma_bar_vals,phi_vals);
% surf(sigma_mesh,phi_mesh,squeeze(E_bend_bar(1,:,:)));
% zlabel('$E$')
% % surf(sigma_mesh,phi_mesh,delA_bar);
% % zlabel('$\Delta A$')
% % surf(sigma_mesh,phi_mesh,rend_bar);
% % zlabel('$r_\mathrm{end}$')
% xlabel('$\bar{\Sigma}$')
% ylabel('$\phi$')
% axes1.XScale = 'log';
% axes1.ZScale = 'log';
% % zlim([1e-5,max(E_bend_bar,[],'all')])
% % ylim([3,pi])

%% test points
% first make some test points
sigma_test = logspace(-5,8,12);
% phi_vals = linspace(0.01,pi-0.1,100);
%phi_vals = deg2rad(linspace(1,179,10));
phi_test = smoothstep(linspace(0.03,1-0.03,7),3)*pi;
d_test = logspace(log10(3.7442), 5,3);

ld = length(d_test);
lp = length(phi_test);
ls = length(sigma_test);

delA_test = nan(ld,lp,ls);
E_bend_test = nan(size(delA_test));
rend_test = nan(size(delA_test));
phi_end_test = nan(size(delA_test));

tic
% ticBytes(gcp);
% parfor (dd = 1:ld, 8)
for dd = 1:ld
    d = d_test(dd)*R;
    for ss = 1:ls
        sigma_bar = sigma_test(ss);
        if mod(ss,1)==0
            fprintf('%i \n', ss);
        end
        for pp = 1:lp
            phi = phi_test(pp);
%             fprintf('%0.4g \n', phi);
            Sigma = sigma_bar*kappa/R^2;
            out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, false);
            r = out.y(2,:);
            rend = r(end);
            rphi = r(1);
            A = out.y(8,end);
            phi_end_test(dd,pp,ss) = out.y(1,end);        
            delA_test(dd,pp,ss) = (2*pi*A - pi*(rend.^2-rphi.^2))/R^2;
            E_bend_test(dd,pp,ss) = out.y(7,end);
            rend_test(dd,pp,ss) = rend/R;
        end
    end
end
% tocBytes(gcp)
toc

%% interpolate at test points
[sigma_mesh,phi_mesh, d_mesh] = meshgrid(sigma_test,phi_test,d_test);
[smi,pmi,dmi] = meshgrid(sigma_bar_vals,phi_vals,d_bar_vals);
phi_interp = permute(interp3(smi,pmi,dmi,permute(phi_end,[2,3,1]),sigma_mesh,phi_mesh,d_mesh),[3,1,2]);
delA_interp = permute(interp3(smi,pmi,dmi,permute(delA_bar,[2,3,1]),sigma_mesh,phi_mesh,d_mesh),[3,1,2]);
rend_interp = permute(interp3(smi,pmi,dmi,permute(rend_bar,[2,3,1]),sigma_mesh,phi_mesh,d_mesh),[3,1,2]);
E_bend_interp = permute(interp3(smi,pmi,dmi,permute(E_bend_bar,[2,3,1]),sigma_mesh,phi_mesh,d_mesh),[3,1,2]);

% phi_interp_perm = permute(phi_interp,[3,1,2]);
rel_E = (E_bend_interp-E_bend_test)./E_bend_test;

squeeze(E_bend_test(:,:,1))
squeeze(E_bend_interp(:,:,1))

