%% combine variables back together
num_tasks = 48;
results = [];
for ii=1:num_tasks-1
    load_name = sprintf('data/task%i_results.mat', ii );
    data = load(load_name);
    results = cat(2,results,data.results_self_mat);
end

delA_bar = squeeze(results(1,:,:,:));
E_bend_bar = squeeze(results(2,:,:,:));
rend_bar = squeeze(results(3,:,:,:));
psi_end = squeeze(results(4,:,:,:));

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
phi_test = smoothstep(linspace(0.03,1-0.03,8),3)*pi;
m_test = [-1,0,1];

lm = length(m_test);
lp = length(phi_test);
ls = length(sigma_test);

delA_test = nan(lm,lp,ls);
E_bend_test = nan(size(delA_test));
rend_test = nan(size(delA_test));
psi_end_test = nan(size(delA_test));

tic
% ticBytes(gcp);
% parfor (dd = 1:ld, 8)
for mm = 1:lm
    m = m_test(mm)*R;
    for ss = 1:ls
        sigma_bar = sigma_test(ss);
        if mod(ss,1)==0
            fprintf('%i \n', ss);
        end
        for pp = 1:lp
            phi = phi_test(pp);
%             fprintf('%0.4g \n', phi);
            Sigma = sigma_bar*kappa/R^2;
            out = Canalejo_2021_m_func(R, phi, Sigma, kappa, m, d, false);
            r = out.y(2,:);
            rend = r(end);
            rphi = r(1);
            delA_test(mm,pp,ss) = out.y(6,end)/R^2;
            E_bend_test(mm,pp,ss) = out.y(7,end)/kappa ...
                - 2*pi*(rend^2-rphi^2)*m^2;
            rend_test(mm,pp,ss) = rend/R;
            psi_end_test(mm,pp,ss) = out.y(1,end);
        end
    end
end
% tocBytes(gcp)
toc

%% interpolate at test points
[sigma_mesh,phi_mesh, m_mesh] = meshgrid(sigma_test,phi_test,m_test);
[smi,pmi,mmi] = meshgrid(sigma_bar_vals,phi_vals,m_bar_vals);
psi_interp = permute(interp3(smi,pmi,mmi,permute(psi_end,[2,3,1]),...
    sigma_mesh,phi_mesh,m_mesh),[3,1,2]);
delA_interp = permute(interp3(smi,pmi,mmi,permute(delA_bar,[2,3,1]), ...
    sigma_mesh,phi_mesh,m_mesh),[3,1,2]);
rend_interp = permute(interp3(smi,pmi,mmi,permute(rend_bar,[2,3,1]), ...
    sigma_mesh,phi_mesh,m_mesh),[3,1,2]);
E_bend_interp = permute(interp3(smi,pmi,mmi,permute(E_bend_bar,[2,3,1]), ...
    sigma_mesh,phi_mesh,m_mesh),[3,1,2]);

% phi_interp_perm = permute(phi_interp,[3,1,2]);
rel_E = (E_bend_interp-E_bend_test)./E_bend_test;

squeeze(E_bend_test(:,:,1))
squeeze(E_bend_interp(:,:,1))