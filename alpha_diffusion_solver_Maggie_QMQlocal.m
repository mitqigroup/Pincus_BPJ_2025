% clear variables
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');

% lookup_table = load("lookup_table.mat")
lookup_table = load("../Scripts/Membranes_Blood/Nonlinear_solutions/gridded_interp.mat");
% lookup_table = load("gridded_interp.mat");

%addpath("/home/gridsan/ipincus/membranes/codes_scripts/Membrane_scripts");
addpath("../Scripts/Membranes_Blood/Linear_solution");
addpath("../Scripts/Membranes_Blood" + ...
    "/Nonlinear_solutions");

MyTaskID = 48;
NumberOfTasks = 10;

% check that the environment variables have been read in correctly
if ~(exist('MyTaskID', 'var')&&exist('NumberOfTasks', 'var'))
    error('Environment variables not set correctly')
end

% taskIDs count from zero, alter this here
%MyTaskID = MyTaskID + 1;
MyTaskID = 25;
% make a folder and go into it
dir_name = sprintf('task%s_kappa20', string(MyTaskID));
mkdir(dir_name);
cd(dir_name);

ep_vals = -linspace(0.07,0.2,8)/2
alpha_i0_vals = linspace(0.5,2,6)*1e-5;

[m,n] = ndgrid(ep_vals, alpha_i0_vals);
all_vals = [m(:),n(:)];

my_vals = all_vals(MyTaskID, :);

epsilon = my_vals(1);
alpha_i0 = my_vals(2);

% firstly set up parameters, all in units of kBT and nm, so we have that
% approx J/m^2 \equiv 250 kBT/nm^2;
% kD = 0.3*250;
kD = 1*250;
kappa = 20;
% epsilon = -6e-4*250;
n0 = 1;
% R = 25;
% d = 500;
d = 2500; %should have been a variable
% alpha_i0 = -3e-2;
% alpha_i0 = 1.2e-5;
N = 1e4;
A_cell = 800e6;
N_sites = A_cell/d^2;
dt = 0.05;
tau = 3.4;

% R_vals = [12,15,18,20,23,25,28,30,40,50];
%R_vals = [15,16,17,18,20,22,24,25,30,34,43,50];
R_vals = [20,22,24,25,30,34,43,50];
%R_vals = [34];
% R_vals = 25;

%%
for rv = 1:length(R_vals)
    R = R_vals(rv)
    ii = 1;
    clear alpha_i alpha_A alpha_B phi E SA SB Si S0 S alpha
    alpha_i(ii) = alpha_i0;
    alpha_A(ii) = 0.01;
    alpha_B(ii) = 0.01;
    phi(ii) = pi/6;
    const = [epsilon, n0, d, R, kD, kappa, alpha_i(ii), N];
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        get_nonlinear_minimum_lookup_phi(const, [alpha_A(ii), ...
        alpha_B(ii), phi(ii)], ...
        lookup_table);
    alpha_A(ii) = out(1);
    alpha_B(ii) = out(2);
    phi(ii) = out(3);
    [E(ii), SA(ii), SB(ii)] = get_energy(alpha_A(ii), alpha_B(ii), ...
        [const, phi(ii)]);
    Si(ii) = d^2;
    S0(ii) = Si(ii)/(1+alpha_i(ii));
    S(ii) = SA(ii)+SB(ii);
    alpha(ii) = S(ii)/S0(ii)-1;
    tic
    while true
        ii = ii+1;
        if mod(ii,50)==0
            disp([ii phi(ii-1)]);
        end
        d_alpha = -dt/tau*(alpha(ii-1)-alpha_i0);
        alpha(ii) = alpha(ii-1) + d_alpha;
        dS0 = S(ii-1)/(1+alpha(ii))-S(ii-1)/(1+alpha(ii-1));
        S0(ii) = S0(ii-1) + dS0;
    %     Si(ii) = Si(ii-1)+dS0;
        Si(ii) = d^2;
    %    d = sqrt(Si(ii));
        alpha_i(ii) = Si(ii)/S0(ii)-1;
    
        % now solve for the new values
        const = [epsilon, n0, d, R, kD, kappa, alpha_i(ii), N];
    %     [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    %         get_nonlinear_minimum_lookup_phi(const, [alpha_A(ii-1)*1.001, ...
    %         alpha_B(ii-1)*1.001, phi(ii-1)/1.001], lookup_table);
        [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
            get_nonlinear_minimum_lookup_phi(const, [alpha_A(ii-1), ...
            alpha_B(ii-1), phi(ii-1)], lookup_table);
        alpha_A(ii) = out(1);
        alpha_B(ii) = out(2);
        phi(ii) = out(3);
        [E(ii), SA(ii), SB(ii)] = get_energy(alpha_A(ii), alpha_B(ii), ...
            [const, phi(ii)]);
        S(ii) = SA(ii)+SB(ii);
        alpha(ii) = S(ii)/S0(ii)-1;
        if phi(ii)>0.95*pi%2.7
            disp("break due to complete wrapping")
            break
        end
        if phi(ii)<1e-3||ii*dt>1000
            break
        end
    end
    toc
    
    if phi(ii)<1e-3
        % here it won't ever wrap. continue to next.
        tw1(rv) = nan;
        tw2(rv) = -log(0.1)*tau;
        continue
    end
    
    %%
%     figure();
%     hold on
%     plot((1:ii)*dt, phi/pi, 'DisplayName',sprintf('$\\Delta t = %0.3g$', dt));
%     % plot((1:ii)*dt, alpha_B, 'DisplayName',sprintf('$\\Delta t = %0.3g$', dt));
%     % plot((1:ii-1)*dt, phi/pi, 'DisplayName',sprintf('$\\Delta t = %0.3g$', dt));
%     % plot((1:ii-1)*dt, alpha_B, 'DisplayName',sprintf('$\\Delta t = %0.3g$', dt));
%     xlabel('time (s)');
%     ylabel('$\phi/\pi$');
%     % axes1 = gca;
%     % axes1.XScale = 'log';
%     legend
%     
    
    %%
    dS_total = S(end)-S(1);
    dS0_total = S0(end)-S0(1);
    dSi_total = Si(end)-Si(1);
    area_consumed = 4*pi*R^2;
    S_particle = area_consumed;
    
    dttau = dt/tau;
    Sdt = S(end);
    S0dt = S0(end);
    adt = alpha(end);
    aidt = alpha_i(end);
    
    Sitw = Si(end);
    S0tw = Sitw/(1+alpha_B(end));
    ai0 = alpha_i(1);
    Stw = Sitw;
    
    t1 = (1:ii)*dt;
    tw1(rv) = ii*dt;
    tw2_lin = (Sitw/(1+ai0)-S0tw)/(S0tw^2/Sitw/tau*(Sitw/S0tw-(1+ai0)));
    % tw2 = -tau*log((Sitw/(0.1*S0tw-0.9*Sitw)-(1+ai0))/(Sitw/S0tw-(1+ai0)));
    tw2(rv) = -tau*log((-Sitw/(0.1*S0tw+0.9*S0(1))+(1+ai0))/(-Sitw/S0tw+(1+ai0)));
    rate_w2 = (S0tw^2/Sitw/tau*(Sitw/S0tw-(1+ai0)));
    
%     figure();
%     hold on
%     plot(t1, S0-Si(1), 'r-', ...
%         'DisplayName','Wrapping, $t_{w,1}$');
%     
    t2 = linspace(0,10,1000);
    S02 = Sitw./((Sitw/S0tw-(1+ai0))*exp(-t2/tau)+1+ai0);
    S02_lin = (t2)*rate_w2+S02(1);
%     plot(t2+tw1(rv), S02-Sitw, 'b-', ...
%         'DisplayName','Replenishment, $t_{w,2}$, Full solution');
%     plot(t2+tw1(rv), S02_lin-Sitw, 'b:', ...
%         'DisplayName','Replenishment, $t_{w,2}$, Linear');
%     plot([t1,t2+tw1(rv)], ones(size([t1,t2]))*(S0(1)-Sitw), 'k--',...
%         'DisplayName','Initial $S_0$')
%     
%     xlabel('time (s)')
%     ylabel('$S_0-S_i$ (nm$^2$)')
    
    legend
    
    % alpha plot
%     figure();
%     hold on
%     plot(t1, alpha, 'r-', ...
%         'DisplayName','Wrapping, $t_{w,1}$');
    
    alpha2 = Sitw./S02-1;
    alpha2_lin = Sitw./S02_lin-1;
%     plot(t2+tw1(rv), alpha2, 'b-', ...
%         'DisplayName','Replenishment, $t_{w,2}$, Full solution');
%     plot(t2+tw1(rv), alpha2_lin, 'b:', ...
%         'DisplayName','Replenishment, $t_{w,2}$, Linear');
    
%     plot([t1,t2+tw1(rv)], ones(size([t1,t2]))*ai0, 'k--',...
%         'DisplayName','Initial $S_0$')
    
%     xlabel('time (s)')
%     ylabel('$\alpha$')
%     
%     legend
    
    % we should have that:
    dS0_pB = S0(end)-(S(end)-S_particle)/(1+alpha_B(end));
    dS0_pA = S_particle/(1+alpha_A(end));
    % are equal! This indeed seems to be the case to a very good approximation.
end

save

% %% overall extraction
% figure();
% hold on
% rate_exp = [0.5, 1.5, 2.5, 2, 1.25]*1e3;
% size_p = R_vals*2;
% rate = 1./(tw1+tw2);
% % uptake = rate/max(rate)*max(rate_exp)*1.3;
% uptake = rate*2*3600*N_sites;
% plot(size_p, uptake, 'ro-', 'color', 'b', 'MarkerFaceColor','b',...
%     'Displayname', 'Model','MarkerSize', 12);
% size_p_exp = [14,30,50,74,100];
% % rate_exp = [3, 4.5, 6, 4, 1.8];
% % rate_exp = rate_exp/max(rate_exp)s;
% plot(size_p_exp, rate_exp, 'rs:', 'color', 'r', 'MarkerFaceColor','r',...
%     'DisplayName', 'Experiments','MarkerSize', 12);
% xlabel('Particle Size');
% ylabel('Particles/cell at 2hrs');
% legend
% % ylim([0,1.1])
% 
% % a1 = annotation('textbox', 'String',...
% %     [sprintf('$R = %0.3g$ nm \n', R),...
% %     sprintf('$2 \\sqrt{2 \\kappa/\\epsilon n_0} = %0.3g $ nm \n', ...
% % 2*sqrt(2*kappa/-epsilon)),...
% %     sprintf('$\\kappa = %0.3g$ $k_\\mathrm{B}T$ \n', kappa),...
% %     sprintf('$d = %0.3g$ nm \n', d),...
% %     sprintf('$\\tau = %0.3g$ s \n', 1/D),...
% %     sprintf('$k_D = %0.3g$ $k_\\mathrm{B}T$/nm$^2$', kD)]);
% a1 = annotation('textbox', 'String',...
%     [sprintf('$2 \\sqrt{2 \\kappa/\\epsilon n_0} = %0.3g $ nm \n', ...
%     2*sqrt(2*kappa/-epsilon)),...
%     sprintf('$\\kappa = %0.3g$ $k_\\mathrm{B}T$ \n', kappa),...
%     sprintf('$d = %0.3g$ nm \n', d),...
%     sprintf('$\\tau = %0.3g$ s \n', tau),...
%     sprintf('$k_D = %0.3g$ $k_\\mathrm{B}T$/nm$^2$', kD)]);
% a1.Position = [0.161266666666667,0.456403269754769, ...
%     0.379246153846154,0.218397820163489];
% 
% %%
% opts = detectImportOptions("Plots\Lipowsky_uptake.csv");
% data = readmatrix("Plots\Lipowsky_uptake.csv", opts);
% 
% p = 1-1e-1;
% X = data(:,1);
% Y = data(:,2)*(1/1.5e-3)*N_sites;
% % Y = data(:,2);
% pp = csaps(X,Y,p);
% % pp = csaps(X,Y);
% fnplt(pp)
% % plot(data(:,1), data(:,2)*(1/1.5e-3))
% % 
% opts = detectImportOptions("Plots\Gao_data.csv");
% data = readmatrix("Plots\Gao_data.csv", opts);
% 
% % p = 0.2;
% X = data(:,1)*(sqrt(kappa/5e-3))*2;
% Y = 2*3600./(data(:,2)*(kappa/(5e3*1e-2)))*N_sites;
% % plot(X,Y, 'rx');
% 
% pp = csaps(X,Y,p);
% % pp = csaps(X,Y);
% fnplt(pp)
% 
% axes1 = gca;
% axes1.YScale = 'log';

%% functions

function [value,isterminal,direction] = myEventsFcn(t,y, cutoff)
    value = y(1)-cutoff;
    isterminal = 1;
    direction = 0;
end

function [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_nonlinear_minimum_lookup_phi(constants, inputs, lookup_table)

    epsilon = constants(1);
    n0 = constants(2);
    d = constants(3);
    R = constants(4);
    kD = constants(5);
    kappa = constants(6);
    alpha_i = constants(7);
    N = constants(8);
    
    alpha_A_init = inputs(1);
    alpha_B_init = inputs(2);
    phi_init = inputs(3);

    phi = phi_init;
    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = 0;
    S_i = d^2;

    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, ...
        'algorithm', 'sqp',...
        'OptimalityTolerance', 1e-14, 'ConstraintTolerance', 1e-14, ...
        'StepTolerance', 1e-14, 'Display', 'off', 'ScaleProblem', true,...
        'TolConSQP', 1e-12, 'FiniteDifferenceType', 'central',...
        'FiniteDifferenceStepSize', sqrt(eps)/1e3);
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init, phi_init],...
        [],[],[],[],[-0.1,alpha_i*0.99, phi_init/2],[5, 5, pi-1e-5], ...
        @constraint, options);

    function [delA, E_bend_free, rend] = get_lookup(Sigma)
        E_bend_free = lookup_table.F_Ebend(Sigma*R^2/kappa,phi,d/R)*pi*kappa;
        delA = lookup_table.F_delA(Sigma*R^2/kappa,phi,d/R)*R^2;
        rend = lookup_table.F_rend(Sigma*R^2/kappa,phi,d/R)*R;
    end

    function f = objective(y_obj)

        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        phi = y_obj(3);
        
        Sigma = kD*alpha_B;

        [delA, E_bend_free, ~] = get_lookup(Sigma);
        if isnan(E_bend_free)
%         if true
            % if it doesn't work, do old method
            
            out_shape = free_shape_nonlinear_free_h(R, d, phi, kappa, ...
                Sigma, false);
            
            S_A = 2*pi*R^2*(1-cos(phi));
            S_B = out_shape.y(8,end)*2*pi + pi*((d/2)^2-out_shape.y(2,end)^2) ...
                + d^2*(1-pi/4);
            
            % stretching, adhesion and bending energy
            f = epsilon*n0*S_A./(1+alpha_A) ...
              + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
              + out_shape.y(7,end)*kappa*pi + 4*pi*kappa*(1-cos(phi));
        else
            % we looked up correctly, get energies etc
            S_A = 2*pi*R^2*(1-cos(phi));
            S_B = delA+(d^2-pi*(R*sin(phi))^2);
            
            % stretching, adhesion and bending energy
            f = epsilon*n0*S_A./(1+alpha_A) ...
              + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
              + E_bend_free + 4*pi*kappa*(1-cos(phi));
        end
%         fprintf('hello \n');
%         out_shape.y(7,end)*kappa*pi
%         E_bend_free
%         delA+(d^2-pi*(R*sin(phi))^2)-d^2
%         S_B-d^2

    end

    function [c,ceq]= constraint(~)
    
        c(1) = 0;
        ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
%         c(2) = 0;
%         ceq(2) = alpha_B-alpha_A;
    
    end

end

function [E,S_A,S_B] = get_energy(alpha_A, alpha_B, constants)

    epsilon = constants(1);
    n0 = constants(2);
    d = constants(3);
    R = constants(4);
    kD = constants(5);
    kappa = constants(6);
    alpha_i = constants(7);
    N = constants(8);
    phi = constants(9);

    % get the shape of the free region
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
        
    out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, false);

    solution = deval(out, linspace(0,out.x(end), 1000));
    r_nonlin = solution(2,:);
    h_nonlin = solution(3,:)-solution(3,end);

    if out.x(end)<d/2
        r_nonlin = [solution(2,1:end-1),d/2];
        h_nonlin = [solution(3,1:end-1),solution(3,end)]-solution(3,end);
    elseif out.x(end)>d/2
        r_nonlin(r_nonlin>d/2) = d/2;
        h_nonlin(r_nonlin>d/2) = 0;
    end
    
    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = out.y(8,end)*2*pi + pi*((d/2)^2-out.y(2,end)^2) + d^2*(1-pi/4);

    E_adhesion = epsilon*n0*S_A./(1+alpha_A);
    E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
    E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
    E_bend_B = out.y(7,end)*kappa*pi;
    E_bend_A = 4*pi*kappa*(1-cos(phi));
    E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;

    E = E-kD/2*(alpha_i.^2*d^2./(1+alpha_i));

end

function [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_nonlinear_minimum_lookup(constants, inputs, lookup_table)

    epsilon = constants(1);
    n0 = constants(2);
    d = constants(3);
    R = constants(4);
    kD = constants(5);
    kappa = constants(6);
    alpha_i = constants(7);
    N = constants(8);
    phi = constants(9);
    
    alpha_A_init = inputs(1);
    alpha_B_init = inputs(2);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = 0;
    S_i = d^2;

    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;

    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, ...
        'algorithm', 'sqp',...
        'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, ...
        'StepTolerance', 1e-12, 'Display', 'off');
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init],...
        [],[],[],[],[-0.1,alpha_i*0.99],[5, 5], ...
        @constraint, options);

    function [delA, E_bend_free, rend] = get_lookup(Sigma)
        E_bend_free = lookup_table.F_Ebend(Sigma*R^2/kappa,phi,d/R)*pi*kappa;
        delA = lookup_table.F_delA(Sigma*R^2/kappa,phi,d/R)*R^2;
        rend = lookup_table.F_rend(Sigma*R^2/kappa,phi,d/R)*R;
    end

    function f = objective(y_obj)

        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        
        Sigma = kD*alpha_B;

        [delA, E_bend_free, ~] = get_lookup(Sigma);
        if isnan(E_bend_free)
            % if it doesn't work, do old method
            
            out_shape = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, false);
            
            S_A = 2*pi*R^2*(1-cos(phi));
            S_B = out_shape.y(8,end)*2*pi + pi*((d/2)^2-out_shape.y(2,end)^2) ...
                + d^2*(1-pi/4);
            
            % stretching, adhesion and bending energy
            f = epsilon*n0*S_A./(1+alpha_A) ...
              + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
              + out_shape.y(7,end)*kappa*pi + 4*pi*kappa*(1-cos(phi));
        else
            % we looked up correctly, get energies etc
            S_A = 2*pi*R^2*(1-cos(phi));
            S_B = delA+(d^2-pi*(R*sin(phi))^2);
            
            % stretching, adhesion and bending energy
            f = epsilon*n0*S_A./(1+alpha_A) ...
              + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
              + E_bend_free + 4*pi*kappa*(1-cos(phi));
        end

    end

    function [c,ceq]= constraint(~)
    
        c(1) = 0;
        ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
        c(2) = 0;
        ceq(2) = alpha_B-alpha_A;
    
    end

end


% % useful set of parameters
% % firstly set up parameters, all in units of kBT and nm, so we have that
% % approx J/m^2 \equiv 250 kBT/nm^2;
% kD = 0.3*250;
% kappa = 10;
% epsilon = -0.5e-3*250;
% n0 = 1;
% % R = 15;
% d = 500;
% % alpha_i0 = -3e-2;
% alpha_i0 = 1e-4;
% N = 1e4;
% 
% % R_vals = [12,15,18,20,23,25,28,30,40,50];
% R_vals = [13,14,16,18,20,22,23,25,30,34,41,50];