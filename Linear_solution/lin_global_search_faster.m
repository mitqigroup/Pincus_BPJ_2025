clear variables

% constants
R = 0.3;                  % um
sigma = 0.01;          % surface fraction
d = sqrt(R^2/sigma);    % um
% phi = pi/12;
kD  = 300/10^12*1e9;    % picoJ/um^2
zeta = 0.02;            % dimensionless
epsilon = -zeta*kD;     % picoJ/um^2
% epsilon = -1;
n0 = 1;                 % fraction
kappa = 1e-19*1e12;     % picoJ
% kappa = 1e-17*1e12;     % picoJ
alpha_i = 0.01;        % fraction
zeta = epsilon*n0/kD;   % dimensionless
N = 3e3;                % number of points in quadrature
plotfigs = 0;

% param_1_vals = logspace(-21,-15, 60)*1e12;
% param_1_vals = logspace(-4,0, 40);
% param_1_vals = logspace(-3,-1, 1);
param_1_vals = 0.02;

for ii = 1:length(param_1_vals)
%     kappa = param_1_vals(ii);
%     zeta = param_1_vals(ii);
%     epsilon = -zeta*kD;     % picoJ/um^2
%     kD = param_1_vals(ii);
%     kappa = 3.333e-7*kD;
    sigma = param_1_vals(ii);
    d = sqrt(R^2/sigma);    % um

    % solve microplastics for initial stretch of A and B
    options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e4);
    const = [zeta, alpha_i, d, R];
    [out,fval,~,~,lam_vals,grad,hessian] = ...
        fmincon(@(y) free_phi(y, const),[alpha_i, alpha_i, 0.3], [],[],[],[],...
        [-1,-1,0],[Inf, Inf, pi/2], ...
        @(y) lipid_con_phi(y,const), options);
    alpha_A_init = out(1);
    alpha_B_init = out(2);
    phi_init = out(3);
    h_phi_init = 0;
    S_A_init = 2*pi*R^2*(1-cos(phi_init));
    S_B_init = d^2 - pi*R^2*sin(phi_init)^2;

    E_adhesion_init = epsilon*n0*S_A_init ./(1+alpha_A_init );
    E_stretch_A_init  = kD/2*(alpha_A_init .^2*S_A_init ./(1+alpha_A_init ));
    E_stretch_B_init  = kD/2*(alpha_B_init .^2*S_B_init ./(1+alpha_B_init ));
    E_total_init = E_adhesion_init+E_stretch_A_init+E_stretch_B_init;
    
    % use initial stretch to solve for minimum attachment

    const = [epsilon, n0, d, R, kD, kappa, alpha_i, N];
    opts = optimoptions(@fmincon,'Algorithm',"sqp");
    rng default % For reproducibility
    problem = createOptimProblem("fmincon",...
        'x0', [alpha_A_init, alpha_B_init, phi_init, h_phi_init],...
        'objective', @(y) stretch_bend_min(y, const),...
        'lb', [-0.1,-0.1,0,-2*d],...
        'ub', [0.1, 0.1, pi/2, 2*d],...
        'nonlcon', @(y) lipid_con_bend(y,const),...
        'options', opts);
    % gs = GlobalSearch;
    st = alpha_i*2;
%     ptmatrix = [alpha_A_init, alpha_B_init, phi_init, h_phi_init;...
%         -st, -st, 5, 0;-st, st, 5, 0;st, -st, 5, 0;st, st, 5, 0;...
%         -st, -st, 5, R;-st, -st, 5, -R;st, st, 5, d/2;st, st, 5, -d/2;...
%         -st, -st, 85, 0;-st, st, 85, 0;st, -st, 85, 0;st, st, 85, 0;...
%         -st, -st, 85, R;-st, -st, 85, -R;st, st, 85, d/2;st, st, 85, -d/2;...
%         -st, -st, 45, 0;-st, st, 45, 0;st, -st, 45, 0;st, st, 45, 0;...
%         -st, -st, 45, R;-st, -st, 45, -R;st, st, 45, d/2;st, st, 45, -d/2];
    alpha_A_vals_inp = [-0.1, -alpha_i, -1e-3, 0, 1e-3, alpha_i, 0.1];
    alpha_B_vals_inp = [-0.1, -alpha_i, -1e-3, 0, 1e-3, alpha_i, 0.1];
    phi_vals_inp = deg2rad([0, 5, 20, 45, 60, 85, 90]);
    h_phi_vals_inp = [-2*d,-d/2, -R, 0, R, d/2,2*d];
    tic
    for aa = 1:length(alpha_A_vals_inp)
%         for bb = 1:length(alpha_B_vals_inp) 
            for pp = 1:length(phi_vals_inp)
                for hh = 1:length(h_phi_vals_inp) 
                    const_con = [epsilon, n0, d, R, kD, kappa, alpha_i, N,...
                        alpha_A_vals_inp(aa), phi_vals_inp(pp), h_phi_vals_inp(hh)];
                    alpha_B_vals_inp(pp,hh,aa) = fzero(@(x) lipid_con_bend_test_A(x, const_con), rand()*0.2-0.1);
                    ptmatrix(aa,pp,hh,1:4) = ...
                        [alpha_A_vals_inp(aa), alpha_B_vals_inp(pp,hh,aa),...
                         phi_vals_inp(pp), h_phi_vals_inp(hh)];
%                     ptmatrix(aa,bb,pp,hh,1:4) = ...
%                         [alpha_A_vals_inp(aa), alpha_A_vals_inp(bb),...
%                          phi_vals_inp(pp), h_phi_vals_inp(hh)];
                end
            end
%         end
    end
    toc
    ptmatrix = reshape(ptmatrix, [numel(ptmatrix)/4, 4]);
    tpoints = CustomStartPointSet(ptmatrix);
    rs = RandomStartPointSet('NumStartPoints',100);
    gs = MultiStart("FunctionTolerance",1e-3, "XTolerance", 1e-3);
    tic
    [out,fval2,exitflag,output,solutions] = run(gs,problem, {tpoints,rs});
    toc
%     [out,fval2,exitflag,output,solutions] = run(gs,problem, {rs});
%     [out,fval2,exitflag,output,solutions] = run(gs,problem, {tpoints});
    
%     figure();
%     hold on
%     plot(arrayfun(@(x)x.Fval,solutions),'k*')
%     xlabel('Solution number')
%     ylabel('Function value')
%     title('Solution Function Values')
    
    % get energies etc
    alpha_A = out(1);
    alpha_B = out(2);
    phi = out(3);
    h_phi = out(4);
    
    % get the shape of the free region
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
    r = linspace(r_phi, d/2,N);
    
    [h,C,S_B, ~, lap_h, hderiv] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi);
    
    S_A = 2*pi*R^2*(1-cos(phi));
    
    E_adhesion = epsilon*n0*S_A./(1+alpha_A);
    E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
    E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
    E_bend_B = kappa/2*2*pi*trapz(r, r.*lap_h.^2);
    E_bend_A = 4*pi*kappa*(1-cos(phi));
    
    E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;
    
    E_all(1,ii) = E;
    E_all(2,ii) = E_adhesion;
    E_all(3,ii) = E_stretch_A;
    E_all(4,ii) = E_stretch_B;
    E_all(5,ii) = E_bend_A;
    E_all(6,ii) = E_bend_B;
    
    alpha_A_vals(ii) = alpha_A;
    alpha_B_vals(ii) = alpha_B;
    phi_vals(ii) = phi;
    h_phi_vals(ii) = h_phi;
    S_A_vals(ii) = S_A;
    S_B_vals(ii) = S_B;
    Sigma_vals(ii) = Sigma;
    
    % plot of shape and nanoparticle
    if plotfigs
        figure('Position',[400,100,700,500]);
        hold on
        axis equal
        xlabel('$r$')
        ylabel('$h$')
        plot(r, h, 'displayname', 'free surface');
        t = linspace(-pi/2,pi/2,1000);
        x = cos(t)*R;
        % y = sin(t)*R+(R*cos(phi)+h(1));
        y = sin(t)*R+R*cos(phi)+h(1);
        plot(x,y, 'displayname', 'microbead')
        % plot(x, (x-sin(phi))*tan(phi)+h(1))
        % legend('location', 'se')
        annotation('textbox', [0.3,0.7,0.4,0.2], 'String',...
            [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
            sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
            sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
            sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
            sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta),...
            sprintf('$\\kappa = %0.2g k_\\mathrm{D} R^2$ \n', kappa/(kD*R^2))], ...
            'interpreter', 'latex','FontSize',16,...
            'FitBoxToText','on','LineStyle','none', 'Color','b')
        annotation('textbox', [0.6,0.7,0.4,0.2], 'String',...
            [sprintf('$\\alpha_A = %0.2g$ \n', alpha_A),...
            sprintf('$\\alpha_B = %0.2g$ \n', alpha_B),...
            sprintf('$\\phi = %0.2g^{\\circ}$ \n', rad2deg(phi)),...
            sprintf('$h_\\phi = %0.2g$ $\\mu$m \n', h_phi)], ...
            'interpreter', 'latex','FontSize',16,...
            'FitBoxToText','on','LineStyle','none', 'Color','k')
        annotation('textbox', [0.55,0.2,0.4,0.2], 'String',...
            [sprintf('$E = %0.3g$ pJ \n', E)], ...
            'interpreter', 'latex','FontSize',16,...
            'FitBoxToText','on','LineStyle','none', 'Color','r')
    end
    
end

% save('kappa_sweep_R300sig01kD300zeta02ai01.mat')

%% plotting solution space
% figure();
% hold on
% ylim([min(cat(1,solutions.Fval))*1.1, 1e-5])
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
% xlabel('Solution number')
% ylabel('Function value')
% % title('Solution Function Values')


% zlim([-1e-3, 1e-5])
for ii = 1:size(solutions,2)
    X0 = solutions(ii).X0{1};
    X = solutions(ii).X;
    alpha_A_val(ii) = X(1);
    alpha_B_val(ii) = X(2);
    phi_val(ii) = X(3);
    h_phi_val(ii) = X(4);
    E_val(ii) = solutions(ii).Fval;
    dist(ii) = vecnorm(out-X);
    E_dist(ii) = E_val(ii) - fval2;
%     plot3(rad2deg(phi_val(ii)), h_phi_val(ii), E_val(ii), 'k.')
end
% figure();
subplot(2,2,1)
hold on
xlim([0,90])
ylim([min(cat(1,solutions.Fval))*1.1, 1e-4])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(rad2deg(phi_val), E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
c_bar = colorbar;
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$\phi$')
ylabel('$E$')

% figure();
subplot(2,2,2)
hold on
xlim([-R, R])
ylim([min(cat(1,solutions.Fval))*1.1, 1e-4])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(h_phi_val, E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
c_bar = colorbar;
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$h_\phi$')
ylabel('$E$')

% figure();
subplot(2,2,3)
hold on
xlim([-0.1, 0.1])
ylim([min(cat(1,solutions.Fval))*1.1, 1e-4])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(alpha_A_val, E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
c_bar = colorbar;
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$\alpha_\mathrm{A}$')
ylabel('$E$')

% figure();
subplot(2,2,4)
hold on
xlim([-0.1, 0.1])
ylim([min(cat(1,solutions.Fval))*1.1, 1e-4])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(alpha_B_val, E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
c_bar = colorbar;
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$\alpha_\mathrm{B}$')
ylabel('$E$')

figure();
hold on
xlim([0,90])
ylim([-R, R])
zlim([min(cat(1,solutions.Fval))*1.1, 1e-4])
scatter3(rad2deg(phi_val), h_phi_val, E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
c_bar = colorbar;
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
set(gca,'ColorScale','log')
xlabel('$\phi$')
ylabel('$h_\phi$')
zlabel('$E$')

% figure();
% hold on
% ylim([0, 1e-3])
% plot(dist, E_dist, '*')
% xlabel('$\Delta X$')
% ylabel('$\Delta E$')

% figure();
% hold on
% xlabel('Solution number')
% ylabel('Multiplicity')
% plot(arrayfun(@(x)length(x.X0),solutions), 'r*')

% %% scatter plot attempt
% figure();
% hold on
% ylim([min(cat(1,solutions.Fval))*1.1, 1e-4])
% xlabel('Solution number')
% ylabel('Function value')
% scatter(1:length(solutions),arrayfun(@(x)x.Fval,solutions),[],...
%     arrayfun(@(x)length(x.X0),solutions), 'filled');
% c_bar = colorbar;
% c_bar.Label.String = 'Multiplicity';
% c_bar.TickLabelInterpreter = 'latex';
% colormap winter
% % set(gca,'ColorScale','log')

%% get surfaces
size_phi = 3;
size_h = 3;
size_A = 3;
phi_test_vals = deg2rad(linspace(eps(1),90-eps(1),size_phi));
h_phi_test_vals = linspace(-2*R, 2*R, size_h);
alpha_A_test_vals = linspace(-0.1,0.1,size_A);
alpha_B_test = zeros(size_phi,size_h,size_A);
E_test = zeros(size_phi,size_h,size_A);
for aa = 1:length(alpha_A_test_vals)
    alpha_A_test = alpha_A_test_vals(aa)
    for pp = 1:length(phi_test_vals)
        phi_test = phi_test_vals(pp);
        for hh = 1:length(h_phi_test_vals)
            h_phi_test = h_phi_test_vals(hh);
            const = [epsilon, n0, d, R, kD, kappa, alpha_i, N,...
                alpha_A_test, phi_test, h_phi_test];
            alpha_B_test(pp,hh,aa) = fzero(@(x) lipid_con_bend_test_A(x, const), rand()*0.2-0.1);
            E_test(pp,hh,aa) = stretch_bend_min([alpha_A_test, alpha_B_test(pp,hh,aa), phi_test, h_phi_test],...
                [epsilon, n0, d, R, kD, kappa, alpha_i, N]);
        end
    end
end

%%
min_phi = 1;
max_phi = 24;
min_h = 1;
max_h = 25;
E_nans = E_test;
E_nans(alpha_B_test<-0.1|alpha_B_test>0.1) = NaN;
figure();
% zlim([-1,1])
hold on
for ii=1:size_A
% surf(rad2deg(phi_test_vals(min_phi:max_phi)),...
%     h_phi_test_vals(min_h:max_h),...
%     squeeze(E_test(min_phi:max_phi,min_h:max_h,1))');
surf(rad2deg(phi_test_vals(min_phi:max_phi)),...
    h_phi_test_vals(min_h:max_h),...
    squeeze(E_nans(min_phi:max_phi,min_h:max_h,ii))',...
    squeeze(alpha_B_test(min_phi:max_phi,min_h:max_h,ii))');
% surf(rad2deg(phi_test_vals(min_phi:max_phi)),...
%     h_phi_test_vals(min_h:max_h),...
%     squeeze(alpha_B_test(min_phi:max_phi,min_h:max_h,ii))');
end

% squeeze(alpha_B_test(min_phi:max_phi,min_h:max_h,1))

[min_val, Q] = min(E_nans, [], 'all');
sE = size(E_nans);
ct = floor(Q/(sE(1)*sE(2)))+1;
bt = floor((Q - (ct-1)*(sE(1)*sE(2)))/sE(1))+1;
at = floor(Q-(ct-1)*(sE(1)*sE(2))-(bt-1)*sE(1));

% phi_min = rad2deg(phi_test_vals(at));
phi_min = phi_test_vals(at);
h_phi_min = h_phi_test_vals(bt);
alpha_A_min = alpha_A_test_vals(ct);
alpha_B_min = alpha_B_test(at,bt,ct);
out_approx_min = [alpha_A_min, alpha_B_min, phi_min, h_phi_min];
% 
plot3(rad2deg(phi_min), h_phi_min, min_val, '.g', 'MarkerSize',30)
plot3(rad2deg(out(3)), out(4), min_val, '.r', 'MarkerSize',30)

c_bar = colorbar;
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = '$\alpha_B$';

set(gca, 'Clim', [0, 0.1])

xlabel('$\phi$')
ylabel('$h_\phi$')
zlabel('$E$')

%% get surfaces alpha B
size_phi = 3;
size_h = 3;
size_B = 3;
phi_test_vals = deg2rad(linspace(eps(1),90-eps(1),size_phi));
h_phi_test_vals = linspace(-2*R, 2*R, size_h);
min_aB = min(alpha_B_test(alpha_B_test>-0.1), [], 'all');
max_aB = max(alpha_B_test(alpha_B_test<0.1), [], 'all');
alpha_B_test_vals = linspace(min_aB,max_aB,size_B);
alpha_A_test = zeros(size_phi,size_h,size_B);
E_test = zeros(size_phi,size_h,size_B);
for bb = 1:length(alpha_B_test_vals)
    alpha_B_test = alpha_B_test_vals(bb)
    for pp = 1:length(phi_test_vals)
        phi_test = phi_test_vals(pp);
        for hh = 1:length(h_phi_test_vals)
            h_phi_test = h_phi_test_vals(hh);
            const = [epsilon, n0, d, R, kD, kappa, alpha_i, N,...
                alpha_B_test, phi_test, h_phi_test];
            alpha_A_test(pp,hh,bb) = fzero(@(x) lipid_con_bend_test_B(x, const), [-.5, 1e10]);
            E_test(pp,hh,bb) = stretch_bend_min([alpha_A_test(pp,hh,bb), alpha_B_test, phi_test, h_phi_test],...
                [epsilon, n0, d, R, kD, kappa, alpha_i, N]);
        end
    end
end

%% solve for alpha_A
testing_vals = linspace(0,10,1000);
alpha_B_test = 0.000011;
phi_test = deg2rad(eps(1)+30);
h_phi_test = -0.0006;
% d = 3;
const = [epsilon, n0, d, R, kD, kappa, alpha_i, N,...
            alpha_B_test, phi_test, h_phi_test];
for bb = 1:length(testing_vals)
%     alpha_B_test = alpha_B_test_vals(bb);
    alpha_A_test = testing_vals(bb);
    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(phi_test));
    func_val(bb) = lipid_con_bend_test_B(alpha_A_test, const);
end

figure();
plot(testing_vals, func_val);

fzero(@(x) lipid_con_bend_test_B(x, const), [-.5, 1e10])

%% solve for alpha_B
testing_vals = linspace(0,0.1,1000);
alpha_A_test = 0.05;
phi_test = deg2rad(eps(1)+20);
h_phi_test = -0.0006;
% d = 3;
const = [epsilon, n0, d, R, kD, kappa, alpha_i, N,...
            alpha_A_test, phi_test, h_phi_test];
for bb = 1:length(testing_vals)
%     alpha_B_test = alpha_B_test_vals(bb);
    alpha_B_test = testing_vals(bb);
    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(phi_test));
    func_val(bb) = lipid_con_bend_test_A(alpha_B_test, const);
end

figure();
plot(testing_vals, func_val);

fzero(@(x) lipid_con_bend_test_A(x, const), [-.5, 1e10])

%%
min_phi = 1;
max_phi = 25;
min_h = 1;
max_h = 25;
E_nans = E_test;
E_nans(alpha_A_test<-0.1|alpha_A_test>0.1) = NaN;
figure();
% zlim([-1,1])
hold on
for ii=1:size_B
% surf(rad2deg(phi_test_vals(min_phi:max_phi)),...
%     h_phi_test_vals(min_h:max_h),...
%     squeeze(E_test(min_phi:max_phi,min_h:max_h,1))');
surf(rad2deg(phi_test_vals(min_phi:max_phi)),...
    h_phi_test_vals(min_h:max_h),...
    squeeze(E_nans(min_phi:max_phi,min_h:max_h,ii))');
end

% squeeze(alpha_B_test(min_phi:max_phi,min_h:max_h,1))

[min_val, Q] = min(E_nans, [], 'all');
sE = size(E_nans);
ct = floor(Q/(sE(1)*sE(2)))+1;
bt = floor((Q - (ct-1)*(sE(1)*sE(2)))/sE(1))+1;
at = floor(Q-(ct-1)*(sE(1)*sE(2))-(bt-1)*sE(1));

% phi_min = rad2deg(phi_test_vals(at));
phi_min = phi_test_vals(at);
h_phi_min = h_phi_test_vals(bt);
alpha_A_min = alpha_A_test_vals(ct);
alpha_B_min = alpha_B_test(at,bt,ct);
out_approx_min = [alpha_A_min, alpha_B_min, phi_min, h_phi_min];

plot3(rad2deg(phi_min), h_phi_min, min_val, '.g', 'MarkerSize',30)
plot3(rad2deg(out(3)), out(4), min_val, '.r', 'MarkerSize',30)

xlabel('$\phi$')
ylabel('$h_\phi$')
zlabel('$E$')

%% plotting
init_vals = param_1_vals(1:20);
on = ones(size(init_vals));

figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
ylabel('E (pJ)')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
for ii=1:6
    plot(param_1_vals, E_all(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
% plot(init_vals, on*E_total_init, 'k-.', 'linewidth', 0.5);
% plot(init_vals, on*E_adhesion_init, 'b-.', 'linewidth', 0.5);
% plot(init_vals, on*E_stretch_B_init, 'r-.', 'linewidth', 0.5);
legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
    '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
    'location', 'best')
annotation('textbox', [0.3,0.6,0.4,0.2], 'String',...
    [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
    sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
    sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
    sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
    sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta)],...
    'interpreter', 'latex','FontSize',16,...
    'FitBoxToText','on','LineStyle','none', 'Color','b')

% positive percentages
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
ylabel('Fraction of Energy')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
E_pos = E_all(1,:)-E_all(2,:);
E_pos_2 = sum(E_all(3:end,:));
E_perc = E_all./E_pos;
for ii=3:6
    plot(param_1_vals, E_perc(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
legend({'$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
    'location', 'best')

% areas
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
yyaxis left
ylabel('Area')
% plot(init_vals, on*S_B_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(init_vals, on*(S_B_init+S_A_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(param_1_vals, S_B_vals, '--', 'displayname', '$S_B$')
plot(param_1_vals, S_A_vals+S_B_vals, '-', 'displayname', '$S_\mathrm{total}$')
yyaxis right
ylabel('Area')
% plot(init_vals, on*S_A_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(param_1_vals, S_A_vals, '--', 'displayname', '$S_A$')
legend

% stretches
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
ylabel('$\alpha$')
plot(param_1_vals, alpha_B_vals, 'displayname', '$\alpha_B$')
plot(param_1_vals, alpha_A_vals, 'displayname', '$\alpha_A$')
% plot(init_vals, on*alpha_A_init, 'r-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(init_vals, on*alpha_B_init, 'b-.', 'linewidth', 0.5, 'HandleVisibility','off');
legend

% phi and h_phi
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
yyaxis left
ylabel('$\phi$')
plot(param_1_vals, rad2deg(phi_vals), 'displayname', '$\phi$')
% plot(init_vals, on*rad2deg(phi_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off')
yyaxis right
ylabel('$h_\phi$')
plot(param_1_vals, h_phi_vals, 'displayname', '$h_\phi$')
legend

%% replotting each function

figure('Position',[400,100,700,500]);
hold on
axis equal
xlabel('$r$')
ylabel('$h$')
for ii = [1,5,10,20,30,40]
    sigma = param_1_vals(ii);
    d = sqrt(R^2/sigma);    % um

    alpha_A = alpha_A_vals(ii);
    alpha_B = alpha_B_vals(ii);
    phi = phi_vals(ii);
    h_phi = h_phi_vals(ii);
    
    % get the shape of the free region
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
    r = linspace(r_phi, d/2,N);
    
    [h,C,S_B, ~, lap_h, hderiv] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi);

    h1 = plot(r, h, 'displayname', sprintf('$\\kappa = %0.2e$', kappa));
    colour = h1.Color;
    t = linspace(-pi/2,-pi/2+phi,1000);
    x = cos(t)*R;
    % y = sin(t)*R+(R*cos(phi)+h(1));
    y = sin(t)*R+R*cos(phi)+h(1);
    plot(x,y, ':', 'HandleVisibility','off', 'color', colour)
end

legend

%% extra functions
function f = stretch_bend_min(y, const)
    % function of free energy to minimise

    epsilon = const(1);
    n0 = const(2);
    d = const(3);
    R = const(4);
    kD = const(5);
    kappa = const(6);
    alpha_i = const(7);
    N = const(8);
    
    alpha_A = y(1);
    alpha_B = y(2);
    phi = y(3);
    h_phi = y(4);

    % problems when phi = 0, so change slightly
    if phi==0
        phi = phi+eps(1);
    end
    
    Sigma = kD*alpha_B;
%     lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
%     r = linspace(r_phi, d/2,N);
    
    [~,S_B, ~,E_bend] = free_shape_linear_no_curve(...
        r_phi, d, phi, kappa, Sigma, h_phi);
    
%     hderiv = C(1)./r+C(3)/lambda*besseli(1,r/lambda)-C(4)/lambda*besselk(1,r/lambda);
%     lap_h = C(3)/lambda^2*besseli(0,r/lambda)+C(4)/lambda^2*besselk(0,r/lambda);
    S_A = 2*pi*R^2*(1-cos(phi));
    
    % stretching, adhesion and bending energy
    f = epsilon*n0*S_A./(1+alpha_A) ...
      + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
      + E_bend + 4*pi*kappa*(1-cos(phi));
% 
%     f = zeta*S_A./(1+alpha_A) ...
%       + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B));

%     f = epsilon*n0*S_A./(1+alpha_A) ...
%       + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B));

end

function [c,ceq] = lipid_con_bend(y, const)

    epsilon = const(1);
    n0 = const(2);
    d = const(3);
    R = const(4);
    kD = const(5);
    kappa = const(6);
    alpha_i = const(7);
    N = const(8);
    
    alpha_A = y(1);
    alpha_B = y(2);
    phi = y(3);
    h_phi = y(4);

    % problems when phi = 0, so change slightly
    if phi==0
        phi = phi+eps(1);
    end
    
    Sigma = kD*alpha_B;
    r_phi = sin(phi)*R;
    
%     r = linspace(r_phi, d/2,N);
    [~,S_B, ~,~] = free_shape_linear_no_curve(...
        r_phi, d, phi, kappa, Sigma, h_phi);
    
    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(phi));

    c(1) = 0;

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end

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


function ceq = lipid_con_bend_test_A(alpha_B, const)

    epsilon = const(1);
    n0 = const(2);
    d = const(3);
    R = const(4);
    kD = const(5);
    kappa = const(6);
    alpha_i = const(7);
    N = const(8);
    alpha_A = const(9);
    phi = const(10);
    h_phi = const(11);

    % problems when phi = 0, so change slightly
    if phi==0
        phi = phi+eps(1);
    end
    
    Sigma = kD*alpha_B;
    r_phi = sin(phi)*R;
    
%     r = linspace(r_phi, d/2,N);
    [~,S_B, ~,~] = free_shape_linear_no_curve(...
        r_phi, d, phi, kappa, Sigma, h_phi);
    
    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(phi));

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end

function ceq = lipid_con_bend_test_B(alpha_A, const)

    epsilon = const(1);
    n0 = const(2);
    d = const(3);
    R = const(4);
    kD = const(5);
    kappa = const(6);
    alpha_i = const(7);
    N = const(8);
    alpha_B = const(9);
    phi = const(10);
    h_phi = const(11);

    % problems when phi = 0, so change slightly
    if phi==0
        phi = phi+eps(1);
    end
    
    Sigma = kD*alpha_B;
    r_phi = sin(phi)*R;
    
%     r = linspace(r_phi, d/2,N);
    [~,S_B, ~,~] = free_shape_linear_no_curve(...
        r_phi, d, phi, kappa, Sigma, h_phi);
    
    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(phi));

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end