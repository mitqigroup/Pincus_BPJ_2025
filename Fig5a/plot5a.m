clear variables
%close all
figure
current_location = './kappa18pt2_kDpt3_tension1_43nm'; %specify file figure 3a
current_location = './kappa18pt2_kDpt3_tension2_43nm'; % SI figures
ds = read_and_unpermute(current_location, 1, 'SortValues', true);
slice = 1:6;
ps = get_param_slices(ds, slice);
newcolors = ["#648FFF","#DC267F","#FFB000","#27E0D8","#1D0C6B","#FE6100","#016D24"];

%% Energy well/barrier and Forces
hold on

opts = detectImportOptions('attachment.csv');
data = readmatrix('attachment.csv', opts);

x = data(:,1);
y = data(:,2);
dy = data(:,4)-y;

e1 = errorbar(x,y,dy, 'ko', 'MarkerFaceColor', 'auto', 'MarkerSize', 10, ...
    'DisplayName','Anselmo 2013');

tic
av = 1; %lower bound Sigma value
rv = 1;
sv = 1:6;
kv = 1;
ev = 1;
pv = 1;
slice = get_slice(rv,sv,kv,ev, pv, av);
outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds, slice, ds.phi_vals,...
    'PlotCurves', false, 'IgnoreNearNaN', false, 'RemoveUnphysical', true);
ps = get_param_slices(ds, slice);
kD = ps.kD_slice(1); %75=0.3/0.004 kT
alpha_i = ps.alpha_i_slice(1);

for ii=1:length(slice)
    R = ps.R_slice(ii);
    d = ps.d_slice(ii);
    phi_loc_min = outputs_analyse.phi_loc_min(ii);
    E_loc_min = outputs_analyse.E_loc_min(ii);
    if ~isempty(phi_loc_min(1))
        phi = phi_loc_min(1);
        phi = phi(1);
        Em = E_loc_min(1);
        E_min(ii) = Em(1) - kD/2*(alpha_i^2*ps.d_slice(ii)^2/(1+alpha_i));
    else
        hphi(ii) = nan;
        E_min(ii) = Eg(ii) - kD/2*(alpha_i^2*ps.d_slice(ii)^2/(1+alpha_i));
        continue
    end
end
toc

%
parameters.V = 57*(1/0.1-1);%corrected value for mice RBC at 10% hematocrit
parameters.A = 91;%corrected value for mice RBC
parameters.Nt = logspace(0,1.5,101);
parameters.R = ps.R_slice(1)/1000;
parameters.sigma_vals = ps.sigma_slice;
sigma_vals=parameters.Nt*(pi*parameters.R^2)/parameters.A;
E_min_vals = interp1(parameters.sigma_vals, E_min, sigma_vals);
expE_min_vals = interp1(parameters.sigma_vals, exp(E_min), sigma_vals);

N_attach = sigma_vals*parameters.A/(pi*parameters.R^2);
xv=sigma_vals./(1-sigma_vals).*expE_min_vals(1);
true_epsilon = ps.epsilon_slice(1)*0.004*10^3
[c index] = min(abs(y(2)-N_attach));
prefactor=(x(2)-N_attach(index))/(xv(index)*parameters.V/(4/3*pi*parameters.R^3));
Nt_vals=xv*prefactor*parameters.V/(4/3*pi*parameters.R^3)+N_attach;
plot(Nt_vals, N_attach, '--', 'DisplayName',...
    'Upper bound, fixed tension','Color', newcolors(av)); 

N_attach = sigma_vals*parameters.A/(pi*parameters.R^2);
xv=sigma_vals./(1-sigma_vals).*expE_min_vals;
true_epsilon = ps.epsilon_slice(1)*0.004*10^3
[c index] = min(abs(y(end)-N_attach));
prefactor=(x(end)-N_attach(index))/(xv(index)*parameters.V/(4/3*pi*parameters.R^3));
Nt_vals=xv*prefactor*parameters.V/(4/3*pi*parameters.R^3)+N_attach;
plot(Nt_vals, N_attach, 'DisplayName',...
    'Upper bound, varying tension',...
    'Color', newcolors(av)); %should be alpha_i * kD
N_attach_1 = N_attach;
Nt_vals1=Nt_vals;

N_attach = sigma_vals*parameters.A/(pi*parameters.R^2);
xv=sigma_vals./(1-sigma_vals).*expE_min_vals(1);
[c index] = min(abs(y(end)-N_attach));
prefactor=(x(end)-N_attach(index))/(xv(index)*parameters.V/(4/3*pi*parameters.R^3));
Nt_vals=xv*prefactor*parameters.V/(4/3*pi*parameters.R^3)+N_attach;
plot(Nt_vals, N_attach, '--', 'DisplayName',...
    'Lower bound, fixed tension','Color', newcolors(av+1)); 

N_attach = sigma_vals*parameters.A/(pi*parameters.R^2);
xv=sigma_vals./(1-sigma_vals).*expE_min_vals;
true_epsilon = ps.epsilon_slice(1)*0.004*10^3
[c index] = min(abs(y(2)-N_attach));
prefactor=(x(2)-N_attach(index))/(xv(index)*parameters.V/(4/3*pi*parameters.R^3));
Nt_vals=xv*prefactor*parameters.V/(4/3*pi*parameters.R^3)+N_attach;
plot(Nt_vals, N_attach, 'DisplayName',...
   'Lower bound, varying tension',...
    'Color', newcolors(av+1)); %should be alpha_i * kD
N_attach_2 = N_attach;
Nt_vals2=Nt_vals;

%% patch
nt_vals_patch = [Nt_vals1, flip(Nt_vals2)];
N_attach_patch = [N_attach_1, flip(N_attach_2)];

p1 = patch(nt_vals_patch,N_attach_patch, [1,1,1]*0.8, 'HandleVisibility', 'off');
uistack(p1, 'bottom');

legend('Position',[0.2,0.71,0.355140575568176,0.155653954461745])

xlim([-1,100.1]);
ylim([-1,32]);
xlabel('Particles/RBC in bulk')
ylabel('Particles/RBC adhered')



function slice = get_slice(rv, sv, kv, ev, pv, av)
% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = 100;                 % nm
sigma_vals = logspace(-5,-0.3, 6);
kD_vals = [0.3]/0.004;                        % kT/nm^2
kappa_vals = [2];%[0.8 0.9];%[0.9 1.0];%%[0.5 1];%[8 10 20 30]; %[2 4 6];%        % kT
epsilon_vals = -kappa_vals/2/(22.5^2);              % kT/nm^2
n0_vals = 1;                                    % fraction
alpha_i_vals = [2e-6]/0.6;                            % fraction
phi_vals = deg2rad(linspace(0.1,179.9,200));

ii = 0;
slice = [];
for rr = 1:length(R_vals)
    for ss = 1:length(sigma_vals)
        for kk = 1:length(kD_vals)
            for ee = 1:length(epsilon_vals)
                for nn = 1:length(n0_vals)
                    for pp = 1:length(kappa_vals)
                        for aa = 1:length(alpha_i_vals)
                            ii = ii+1;
                            if any(rr==rv) && any(ss==sv) && any(kk==kv) && any(ee==ev) && any(rr==rv) && any(pp==pv) && any(aa==av)
                                slice = [slice, ii];
                            end
                        end
                    end
                end
            end
        end
    end
end

end












