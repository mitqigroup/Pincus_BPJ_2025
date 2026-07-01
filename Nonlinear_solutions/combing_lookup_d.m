clear variables

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

%% combine variables back together
tic
num_tasks = 48;
results = [];
for ii=1:num_tasks-1
    load_name = sprintf('data/task%i_results.mat', ii )
    data = load(load_name);
    results = cat(2,results,data.results_self_mat);
end

phi_end = squeeze(results(1,:,:,:));
delA_bar = squeeze(results(2,:,:,:));
E_bend_bar = squeeze(results(3,:,:,:));
rend_bar = squeeze(results(4,:,:,:));
toc

% delA_bar(phi_end>0.2) = nan;
% E_bend_bar(phi_end>0.2) = nan;
% rend_bar(phi_end>0.2) = nan;

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
phi_test = smoothstep(linspace(0.03,1-0.03,20),3)*pi;
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
delA_interp = permute(interp3(smi,pmi,dmi,permute(delA_bar,[2,3,1]),sigma_mesh,phi_mesh,d_mesh, 'linear'),[3,1,2]);
rend_interp = permute(interp3(smi,pmi,dmi,permute(rend_bar,[2,3,1]),sigma_mesh,phi_mesh,d_mesh),[3,1,2]);
E_bend_interp = permute(interp3(smi,pmi,dmi,permute(E_bend_bar,[2,3,1]),sigma_mesh,phi_mesh,d_mesh, 'linear'),[3,1,2]);


z_t = 1-cos(phi_mesh);
other_E = permute(4*z_t+sigma_mesh.*z_t.^2,[3,1,2]);
% E_bend_interp(phi_end_test>0.2) = nan;

% phi_interp_perm = permute(phi_interp,[3,1,2]);
rel_E = (E_bend_interp-E_bend_test)./E_bend_test;
rel_E_full = (E_bend_interp-E_bend_test)./(E_bend_test+other_E);

% squeeze(delA_test(:,:,2))
% squeeze(delA_interp(:,:,2))
rel_A = (delA_interp-delA_test)./delA_test

% squeeze(E_bend_test(:,:,2))
% squeeze(E_bend_interp(:,:,2))

%% test individual points
sig_t = 1e8;
phi_t = 0.0001;
d_t = 600;
Sigma = sig_t*kappa/R^2;
out = free_shape_nonlinear_free_h(R, d_t*R, phi_t, kappa, Sigma, false);
r = out.y(2,:);
rend = r(end);
rphi = r(1);
A = out.y(8,end);
delA_t = (2*pi*A - pi*(rend.^2-rphi.^2))/R^2;
E_bend_t = out.y(7,end);
true_E = E_bend_t
test_E = interp3(smi,pmi,dmi,permute(E_bend_bar,[2,3,1]),sig_t,phi_t,d_t, 'linear')

z = 1-cos(phi_t);
other_E = 4*z+sig_t*z.^2
free_E = (out.y(9,end)-out.y(7,end))

total_rel = (true_E-test_E)/(true_E+other_E)

%% create a gridded interpolant and save it
[sgi,pgi,dgi] = ndgrid(sigma_bar_vals,phi_vals,d_bar_vals);
% [sigma_grid, phi_grid, d_grid] = ndgrid(sigma_test,phi_test,d_test);
F_Ebend = griddedInterpolant(sgi,pgi,dgi,permute(E_bend_bar,[3,2,1]), 'linear');
F_delA = griddedInterpolant(sgi,pgi,dgi,permute(delA_bar,[3,2,1]), 'linear');
F_rend = griddedInterpolant(sgi,pgi,dgi,permute(rend_bar,[3,2,1]), 'linear');
% E_bend_interp_grid = permute(F_Ebend(sigma_grid,phi_grid,d_grid),[3,2,1]);

save("gridded_interp.mat","F_Ebend","F_delA","F_rend");

