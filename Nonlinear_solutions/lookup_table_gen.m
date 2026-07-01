% actually everything depends only on Sigma and phi, but we'll write out
% some example values just to be safe
alpha_i = 1e-4;
R = 0.1;
kappa = 1e-7;
sigma_bar_vals = logspace(1,2,20);
% phi_vals = linspace(0.01,pi-0.1,100);
phi_vals = deg2rad(linspace(0.01,179.9,5));
% phi_vals = smoothstep(linspace(0.001,1-0.001,150),6)*pi;
d = 100000;

delA_bar = nan(length(phi_vals),length(sigma_bar_vals));
E_bend_bar = nan(size(delA_bar));
rend_bar = nan(size(delA_bar));
phi_end = nan(size(delA_bar));

tic
for ss = 1:length(sigma_bar_vals)
    sigma_bar = sigma_bar_vals(ss);
    for pp = 1:length(phi_vals)
        phi = phi_vals(pp);
        Sigma = sigma_bar*kappa/R^2;
        out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, 0);
        r = out.y(2,:);
        rend = r(end);
        rphi = r(1);
        A = out.y(8,end);
        phi_end(pp,ss) = out.y(1,end);        
        delA_bar(pp,ss) = (2*pi*A - pi*(rend.^2-rphi.^2))/R^2;
        E_bend_bar(pp,ss) = out.y(7,end);
        rend_bar(pp,ss) = rend/R;
    end
end
toc


%%
figure();
hold on
axes1 = gca;
[sigma_mesh,phi_mesh] = meshgrid(sigma_bar_vals,phi_vals);
surf(sigma_mesh,phi_mesh,E_bend_bar);
zlabel('$E$')
% surf(sigma_mesh,phi_mesh,delA_bar);
% zlabel('$\Delta A$')
% surf(sigma_mesh,phi_mesh,rend_bar);
% zlabel('$r_\mathrm{end}$')
xlabel('$\bar{\Sigma}$')
ylabel('$\phi$')
axes1.XScale = 'log';
axes1.ZScale = 'log';
zlim([1e-5,max(E_bend_bar,[],'all')])

%% sample at lower accuracy, then do interpolation
sigma_bar_vals_coarse = logspace(-5,8,200);
% z_vals_coarse = linspace(0.0001,2-0.0001,50);
% phi_vals_coarse = flip(deg2rad(linspace(0.01,179.9,150)));
% phi_vals_coarse = acos(1-z_vals_coarse);
phi_vals_coarse = smoothstep(linspace(0.001,1-0.001,150),4)*pi;
d = 100000;

delA_bar_coarse = nan(length(phi_vals_coarse),length(sigma_bar_vals_coarse));
E_bend_bar_coarse = nan(size(delA_bar_coarse));
rend_bar_coarse = nan(size(delA_bar_coarse));
phi_end_coarse = nan(size(delA_bar_coarse));

tic
for ss = 1:length(sigma_bar_vals_coarse)
    sigma_bar = sigma_bar_vals_coarse(ss);
    for pp = 1:length(phi_vals_coarse)
        phi = phi_vals_coarse(pp);
        Sigma = sigma_bar*kappa/R^2;
        out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, 0);
        r = out.y(2,:);
        rend = r(end);
        rphi = r(1);
        A = out.y(8,end);
        phi_end_coarse(pp,ss) = out.y(1,end);
        if phi_end_coarse(pp,ss)<0.01
            delA_bar_coarse(pp,ss) = (2*pi*A - pi*(rend.^2-rphi.^2))/R^2;
            E_bend_bar_coarse(pp,ss) = out.y(7,end);
            rend_bar_coarse(pp,ss) = rend/R;
        end
    end
end
toc

%%
figure();
hold on
axes1 = gca;
[sigma_mesh_coarse,phi_mesh_coarse] = meshgrid(sigma_bar_vals_coarse,phi_vals_coarse);
surf(sigma_mesh_coarse,phi_mesh_coarse,E_bend_bar_coarse);
zlabel('$E$')
% surf(sigma_mesh,phi_mesh,delA_bar);
% zlabel('$\Delta A$')
% surf(sigma_mesh,phi_mesh,rend_bar);
% zlabel('$r_\mathrm{end}$')
xlabel('$\bar{\Sigma}$')
ylabel('$\phi$')
axes1.XScale = 'log';
axes1.ZScale = 'log';
zlim([1e-10,max(E_bend_bar_coarse,[],'all')])

figure();
hold on
axes1 = gca;
[sigma_mesh_coarse,phi_mesh_coarse] = meshgrid(sigma_bar_vals_coarse,phi_vals_coarse);
% surf(sigma_mesh_coarse,phi_mesh_coarse,E_bend_bar_coarse);
zlabel('$E$')
% surf(sigma_mesh,phi_mesh,delA_bar);
% zlabel('$\Delta A$')
% surf(sigma_mesh,phi_mesh,rend_bar);
% zlabel('$r_\mathrm{end}$')
xlabel('$\bar{\Sigma}$')
ylabel('$\phi$')
axes1.XScale = 'log';
axes1.ZScale = 'log';

range_phi = [2:149];
% range_y = []
% E_interp = interp2(sigma_mesh_coarse,phi_mesh_coarse,E_bend_bar_coarse,...
%     sigma_mesh,phi_mesh,'spline');
E_interp = interp2(sigma_mesh_coarse(range_phi,:),phi_mesh_coarse(range_phi,:),E_bend_bar_coarse(range_phi,:),...
    sigma_mesh,phi_mesh,'linear', nan);
% p = reshape(phi_end<0.01,size(sigma_mesh));
% pc = reshape(phi_end_coarse<0.01,size(sigma_mesh_coarse));
% smcn = sigma_mesh_coarse;
% smcn(~pc) = nan;
% pmcn = phi_mesh_coarse;
% pmcn(~pc) = nan;
% emcn = E_bend_bar_coarse;
% emcn(~pc) = nan;
% smn = sigma_mesh;
% smn(~p) = nan;
% pmn = phi_mesh;
% pmn(~p) = nan;
% E_interp = interp2(smcn,pmcn,emcn,smn,pmn,'linear');
surf(sigma_mesh,phi_mesh,E_interp);

lookup_table.sigma_mesh = sigma_mesh_coarse;
lookup_table.phi_mesh = phi_mesh_coarse;
lookup_table.E_bend_bar = E_bend_bar_coarse;
lookup_table.delA_bar_mesh = delA_bar_coarse;
lookup_table.rend_bar = rend_bar_coarse;

rangex = 1:150;
figure();
hold on
axes1 = gca;
% surf(sigma_mesh(rangex,:),phi_mesh(rangex,:),abs((E_interp(rangex,:)-E_bend_bar(rangex,:))./E_bend_bar(rangex,:)));
surf(sigma_mesh,phi_mesh,(E_interp-E_bend_bar)./E_bend_bar);
axes1.XScale = 'log';
axes1.ZScale = 'log';
% zlim([0,0.01])
zlim([-1,1])

% %% interpolate each phi value separately
% 
% E_interp_lin = nan(size(E_interp));
% 
% for pp = 1:length(phi_vals_coarse)
%     E_interp_lin(pp,:) = interp1(sigma_bar_vals_coarse, E_bend_bar_coarse(pp,:), sigma_bar_vals,'cubic');
% end
% 
% figure();
% hold on
% axes1 = gca;
% % [sigma_mesh_coarse,phi_mesh_coarse] = meshgrid(sigma_bar_vals_coarse,phi_vals_coarse);
% % surf(sigma_mesh_coarse,phi_mesh_coarse,E_bend_bar_coarse);
% zlabel('$E$')
% % surf(sigma_mesh,phi_mesh,delA_bar);
% % zlabel('$\Delta A$')
% % surf(sigma_mesh,phi_mesh,rend_bar);
% % zlabel('$r_\mathrm{end}$')
% xlabel('$\bar{\Sigma}$')
% ylabel('$\phi$')
% axes1.XScale = 'log';
% axes1.ZScale = 'log';
% surf(sigma_mesh,phi_mesh,E_interp_lin);
% 
% rangex = 1:150;
% figure();
% hold on
% axes1 = gca;
% surf(sigma_mesh(rangex,:),phi_mesh(rangex,:),abs((E_interp_lin(rangex,:)-E_bend_bar(rangex,:))./E_bend_bar(rangex,:)));
% % surf(sigma_mesh(rangex,:),phi_mesh(rangex,:),abs((E_interp_lin(rangex,:)-E_bend_bar(rangex,:))));
% % surf(sigma_mesh,phi_mesh,(E_interp-E_bend_bar)./E_bend_bar);
% axes1.XScale = 'log';
% axes1.ZScale = 'log';
% % zlim([0,0.01])

