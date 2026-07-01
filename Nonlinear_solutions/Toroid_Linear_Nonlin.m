% clear variables
% close all

% constants
R = 0.2;
sigma = 1e-3;
d = sqrt(pi*R^2/sigma);    % um
% d = 20;    % um
phi = pi/12;
kD  = 300/10^12*1e9;    % picoJ/um^2
kappa = 1e-19*1e12;     % picoJ
kappa_bar = kappa/(kD*R^2);
alpha_i = 0.05;
% zeta = 1e-5;            % dimensionless
% zeta = 1.1*(4*kappa_bar/(2+kappa_bar))          % dimensionless
% zeta = 1.001*(4*kappa_bar*(1+alpha_i))/(2-2*alpha_i-alpha_i^2+4*kappa_bar*(1+alpha_i))
% zeta = 5*(4*kappa_bar*(1+alpha_i))/(2-2*alpha_i-alpha_i^2+4*kappa_bar*(1+alpha_i))
zeta = 0.0000178
% zeta = 0.1*alpha_i^2/2
epsilon = -zeta*kD;     % picoJ/um^2
n0 = 1;                 % fraction
% gamma = 0.1*R;            % um
% epsilon = -kappa/gamma^2;
% kappa = 1e-17*1e12;     % picoJ
% zeta = -epsilon*n0/kD;   % dimensionless
% alpha_i = sqrt(-zeta/2);        % fraction
% alpha_i = 0;
ten_param=-d*epsilon*n0/sqrt(kD*kappa);
psi_dot_init = -1;
N = 3e3;                % number of points in quadrature

savedata = 0;
plotfigs = 0;
filename = '300nm_ai1kap19zeta002.mat';

phi_vals = flip(deg2rad(linspace(0.001,0.2,100)));

% phi_vals = deg2rad(160);

E_all = zeros(6,length(phi_vals));

E_all_toroid = nan(7,length(phi_vals));
alpha_A_vals_toroid = nan(1,length(phi_vals));
alpha_B_vals_toroid = nan(size(alpha_A_vals_toroid));
h_phi_vals_toroid = nan(size(alpha_A_vals_toroid));
S_A_vals_toroid = nan(size(alpha_A_vals_toroid));
S_B_vals_toroid = nan(size(alpha_A_vals_toroid));
Sigma_vals_toroid = nan(size(alpha_A_vals_toroid));
rho_vals = nan(size(alpha_A_vals_toroid));

E_all_linear = nan(7,length(phi_vals));
alpha_A_vals_linear = nan(1,length(phi_vals));
alpha_B_vals_linear = nan(size(alpha_A_vals_toroid));
h_phi_vals_linear = nan(size(alpha_A_vals_toroid));
S_A_vals_linear = nan(size(alpha_A_vals_toroid));
S_B_vals_linear = nan(size(alpha_A_vals_toroid));
Sigma_vals_linear = nan(size(alpha_A_vals_toroid));

E_all_nonlinear = nan(7,length(phi_vals));
alpha_A_vals_nonlinear = nan(1,length(phi_vals));
alpha_B_vals_nonlinear = nan(size(alpha_A_vals_toroid));
h_phi_vals_nonlinear = nan(size(alpha_A_vals_toroid));
S_A_vals_nonlinear = nan(size(alpha_A_vals_toroid));
S_B_vals_nonlinear = nan(size(alpha_A_vals_toroid));
Sigma_vals_nonlinear = nan(size(alpha_A_vals_toroid));

for ii = 1:length(phi_vals)
    phi = phi_vals(ii);
    rad2deg(phi)

    if plotfigs
        figure('Position',[400,100,800,600]);
        hold on
        axis equal
    end
    
    if ii==1
        alpha_B_init = alpha_i;
        alpha_A_init = (1+alpha_B_init)/(1+zeta)-1;
        h_phi_init = 0;
    else
        alpha_A_init = alpha_A;
        alpha_B_init = alpha_B;
        h_phi_init = h_phi;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linear solution, only for phi<85 degrees
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if phi<1.5
        % %% use initial stretch to solve for minimum attachment
        const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi];
            [out,~,~,~,~,~,~] = ...
                get_linear_minimum(const, [alpha_A_init, alpha_B_init, h_phi_init]);
        
        alpha_A = out(1);
        alpha_B = out(2);
        h_phi = out(3);
        
        % get the shape of the free region
        Sigma = kD*alpha_B;
        r_phi = sin(phi)*R;
        r = linspace(r_phi, d/2,N);
        
        [h,C,S_B, ~, hderiv, lap_h] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi);
        
        S_A = 2*pi*R^2*(1-cos(phi));
        
        E_adhesion = epsilon*n0*S_A./(1+alpha_A);
        E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
        E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
        E_bend_B = kappa/2*2*pi*trapz(r, r.*lap_h.^2);
        E_bend_A = 4*pi*kappa*(1-cos(phi));
        
        E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;
        
        E_all_linear(1,ii) = E;
        E_all_linear(2,ii) = E_adhesion;
        E_all_linear(3,ii) = E_stretch_A;
        E_all_linear(4,ii) = E_stretch_B;
        E_all_linear(5,ii) = E_bend_A;
        E_all_linear(6,ii) = E_bend_B;
        
        alpha_A_vals_linear(ii) = alpha_A;
        alpha_B_vals_linear(ii) = alpha_B;
        h_phi_vals_linear(ii) = h_phi;
        S_A_vals_linear(ii) = S_A;
        S_B_vals_linear(ii) = S_B;
        Sigma_vals_linear(ii) = Sigma;
        
        % %% plot of shape and nanoparticle
        if plotfigs
            xlabel('$r$')
            ylabel('$h$')
            plot(r, h, 'r-','displayname', 'free surface');
            t = linspace(-pi/2,pi/2,1000);
            x = cos(t)*R;
            % y = sin(t)*R+(R*cos(phi)+h(1));
            y = sin(t)*R+R*cos(phi)+h(1);
            plot(x,y,'r:', 'displayname', 'microbead')
            % plot(x, (x-sin(phi))*tan(phi)+h(1))
            % legend('location', 'se')
%             annotation('textbox', [0.3,0.7,0.4,0.2], 'String',...
%                 [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
%                 sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
%                 sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
%                 sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
%                 sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta),...
%                 sprintf('$\\kappa = %0.2g$ pJ \n', kappa)], ...
%                 'interpreter', 'latex','FontSize',16,...
%                 'FitBoxToText','on','LineStyle','none', 'Color','b')
%             annotation('textbox', [0.55,0.7,0.4,0.2], 'String',...
%                 [sprintf('$\\alpha_A = %0.2g$ \n', alpha_A),...
%                 sprintf('$\\alpha_B = %0.2g$ \n', alpha_B),...
%                 sprintf('$\\phi = %0.2g^{\\circ}$ \n', rad2deg(phi)),...
%                 sprintf('$h_\\phi = %0.2g$ $\\mu$m \n', h_phi)], ...
%                 'interpreter', 'latex','FontSize',16,...
%                 'FitBoxToText','on','LineStyle','none', 'Color','k')
%             annotation('textbox', [0.55,0.2,0.4,0.2], 'String',...
%                 [sprintf('$E = %0.3g$ pJ \n', E)], ...
%                 'interpreter', 'latex','FontSize',16,...
%                 'FitBoxToText','on','LineStyle','none', 'Color','r')
            % % annotation('textbox', [0.55,0.2,0.4,0.2], 'String',...
            % %     [sprintf('$E = %0.2g$ \n', E),...
            % %     sprintf('$E_\\mathrm{adhesion} = %0.1f\\%%$ \n', E_adhesion/E*100), ...
            % %     sprintf('$E_\\mathrm{stretch,A} = %0.1f\\%%$ \n', E_stretch_A/E*100), ...
            % %     sprintf('$E_\\mathrm{stretch,B} = %0.1f\\%%$ \n', E_stretch_B/E*100), ...
            % %     sprintf('$E_\\mathrm{bend,A} = %0.1f\\%%$ \n', E_bend_A/E*100), ...
            % %     sprintf('$E_\\mathrm{bend,B} = %0.1f\\%%$ \n', E_bend_B/E*100)], ...
            % %     'interpreter', 'latex','FontSize',16,...
            % %     'FitBoxToText','on','LineStyle','none', 'Color','g')
            % 
            % % subplot(1,2,2);
            % axes('Position',[.5 .3 .4 .4])
            % % figure();
            % hold on
            % ylabel('Energy (pJ)')
            % X = categorical({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
            %     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
            %     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'});
            % X = reordercats(X,{'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
            %     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
            %     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'});
            % bar(X, E_all, 'barwidth', 1)
            % set(gca,'xticklabel',{'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
            %     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
            %     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'});
            % bar([-0.5:4.5],E_all, 'barwidth', 1)
            % plot of derivative
            % figure();
            % hold on
            % plot(r,hderiv)
            % plot(r, tan(phi)*ones(size(r)))
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Toroidal solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e5, 'algorithm', 'sqp');
    
%     const = [epsilon, alpha_i, d, R, phi, kappa, kD];
%     [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
%         fmincon(@(y) free_toroid(y, const),[0.001, 0.001, R/20], [],[],[],[],...
%         [-0.1,-0.1,0],[0.1, 0.1,5*d], ...
%         @(y) area_con_toroid(y,const), options);
    

    const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi];
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        get_toroid_minimum(const, [alpha_A_init, alpha_B_init, R/20]);
    
    alpha_A = out(1);
    alpha_B = out(2);
    rho = out(3);
    
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
    
    E = E_adhesion+E_stretch_A+E_stretch_B+E_bend_A+E_bend_B;
    
    E_all_toroid(1,ii) = E;
    E_all_toroid(2,ii) = E_adhesion;
    E_all_toroid(3,ii) = E_stretch_A;
    E_all_toroid(4,ii) = E_stretch_B;
    E_all_toroid(5,ii) = E_bend_A;
    E_all_toroid(6,ii) = E_bend_B;
    
    h_phi = sin(-3*pi/2)*rho+rho*cos(phi)- (sin(-3*pi/2+phi)*rho+rho*cos(phi));
    
    alpha_A_vals_toroid(ii) = alpha_A;
    alpha_B_vals_toroid(ii) = alpha_B;
    h_phi_vals_toroid(ii) = h_phi;
    S_A_vals_toroid(ii) = S_A;
    S_B_vals_toroid(ii) = S_B;
    Sigma_vals_toroid(ii) = alpha_B*kD;
    rho_vals(ii) = rho;
    
    if plotfigs
        hold on

%         t = linspace(-pi/2,-pi/2+phi,1000);
        t = linspace(-pi/2,pi/2,1000);
        x = cos(t)*R;
    
        tt = linspace(-3*pi/2,-3*pi/2+phi,1000);
    %     t = linspace(-pi/2,pi/2,1000);
        xt = cos(tt)*rho;
        % y = sin(t)*R+(R*cos(phi)+h(1));
        yt = sin(tt)*rho+rho*cos(phi);
        plot(xt-xt(end)+cos(-pi/2+phi)*R,yt-yt(1), 'k-', 'HandleVisibility','off')
        
        xf = linspace(xt(1)-xt(end)+cos(-pi/2+phi)*R, d/2,100);
        yf = zeros(size(xf));
        plot(xf, yf, 'k-', 'HandleVisibility','off')


        y = sin(t)*R+R*cos(phi)+(yt(end)-yt(1));
        plot(x,y, 'k:', 'HandleVisibility','off')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     %%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Nonlinear solution
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % %% use initial stretch to solve for minimum attachment
%     % options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e3);
% %     options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp');
% %     psi_dot_init = (rand()-0.5)*3;
% %     const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi, psi_dot_init];
% %     [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
% %         fmincon(@(y) free_nonlinear(y, const),[alpha_A_init, alpha_B_init],...
% %         [],[],[],[],[-0.1,-0.1],[0.1, 0.1], ...
% %         @(y) area_con_nonlinear(y,const), options);
%     
%     const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi];
%     [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
%         get_nonlinear_minimum(const, [alpha_A_init, alpha_B_init]);
%     
%     alpha_A = out(1);
%     alpha_B = out(2);
%     
%     lambda_stretch = lam_vals.eqnonlin;
%     
%     % get the shape of the free region
%     Sigma = kD*alpha_B;
%     lambda = sqrt(kappa/Sigma);
%     r_phi = sin(phi)*R;
%     
%     out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, psi_dot_init);
% 
%     solution = deval(out, linspace(0,out.xe, 1000));
%     r_nonlin = solution(2,:);
%     h_nonlin = solution(3,:)-solution(3,end);
% 
%     if out.x(end)<d/2
%         r_nonlin = [solution(2,1:end-1),d/2];
%         h_nonlin = [solution(3,1:end-1),solution(3,end)]-solution(3,end);
%     elseif out.x(end)>d/2
%         r_nonlin(r_nonlin>d/2) = d/2;
%         h_nonlin(r_nonlin>d/2) = 0;
%     end
%     
%     S_A = 2*pi*R^2*(1-cos(phi));
%     S_B = out.ye(8)*2*pi + pi*((d/2)^2-out.ye(2)^2) + d^2*(1-pi/4);
% 
%     E_adhesion = epsilon*n0*S_A./(1+alpha_A);
%     E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
%     E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
%     E_bend_B = out.ye(7)*kappa*pi;
%     E_bend_A = 4*pi*kappa*(1-cos(phi));
%     E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;
%     
%     E_all_nonlinear(1,ii) = E;
%     E_all_nonlinear(2,ii) = E_adhesion;
%     E_all_nonlinear(3,ii) = E_stretch_A;
%     E_all_nonlinear(4,ii) = E_stretch_B;
%     E_all_nonlinear(5,ii) = E_bend_A;
%     E_all_nonlinear(6,ii) = E_bend_B;
%     
%     alpha_A_vals_nonlinear(ii) = alpha_A;
%     alpha_B_vals_nonlinear(ii) = alpha_B;
%     h_phi_vals_nonlinear(ii) = h_phi;
%     S_A_vals_nonlinear(ii) = S_A;
%     S_B_vals_nonlinear(ii) = S_B;
%     Sigma_vals_nonlinear(ii) = Sigma;
% 
%     % %% plot of shape and nanoparticle
%     % subplot(1,2,1);
%     if plotfigs
% %         figure('Position',[400,100,800,600]);
%         hold on
%         axis equal
%         xlabel('$r$')
%         ylabel('$h$')
%         plot(r_nonlin, h_nonlin, 'b-', 'displayname', 'free surface');
%         t = linspace(-pi/2,pi/2,1000);
%         x = cos(t)*R;
%         % y = sin(t)*R+(R*cos(phi)+h(1));
%         y = sin(t)*R+R*cos(phi)+h_nonlin(1);
%         plot(x,y,'b:', 'displayname', 'microbead')
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

if savedata
    save(filename);
end

%%
% close all
if ~ishandle(1)
    f1 = figure('Position',[400,100,700,500]);
else
    figure(f1);
end
% f1 = figure('Position',[400,100,700,500]);
hold on
% xlim([0,90])
xlabel('$\phi$')
ylabel('$\Delta E$')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
for ii=[1]
%     plot(rad2deg(phi_vals), E_all(ii,:)-E_all(ii,end), ...
%         strcat(colours(ii),lines(ii)))
%     plot(rad2deg(phi_vals), E_all_linear(ii,:), ...
%         strcat(colours(ii),'-'))
    plot(rad2deg(phi_vals), E_all_linear(ii,:)-E_all_linear(ii,end),'-',...
        'displayname', sprintf('$\\zeta = %d$', zeta))
%     plot(rad2deg(phi_vals), E_all_toroid(ii,:), ...
%         strcat(colours(ii),'--'), "HandleVisibility", 'off')
%     plot(rad2deg(phi_vals), E_all_nonlinear(ii,:), ...
%         strcat(colours(ii),':'), "HandleVisibility", 'off')
end
% annotation('textbox', [0.5,0.7,0.4,0.2], 'String',...
%     [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
%     sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
%     sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
%     sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
%     sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta),...
%     sprintf('$\\kappa = %0.2g k_\\mathrm{D} R^2$ \n', kappa/(kD*R^2))], ...
%     'interpreter', 'latex','FontSize',16,...
%     'FitBoxToText','on','LineStyle','none', 'Color','b')
% legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
%     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
%     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
%     'location', 'nw')
% legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
%     '$E_\mathrm{stretch,B}$',...
%     '$E_\mathrm{bend,B}$'}, 'Box','off',...
%     'location', 'nw')
legend

zeta_c = (2*kappa_bar*(1+alpha_i))/(1-alpha_i-1/2*alpha_i^2+2*kappa_bar*(1+alpha_i))

kb = kappa_bar;
ai = alpha_i;

aq = (20-12*ai-4*ai^2+24*kb*(1+ai));
bq = -(-8+8*ai+4*ai^2-16*kb*(1+ai));
cq = (16*kb*(1+ai));

zeta_c2 = (-bq-sqrt(bq^2-4*aq*cq))/(2*aq)
zeta_c3 = (-bq+sqrt(bq^2-4*aq*cq))/(2*aq)

aq = (20+24*kb*(1+ai));
bq = -(-8-16*kb*(1+ai));
cq = (16*kb*(1+ai));
zeta_c5 = (-bq+sqrt(bq^2-4*aq*cq))/(2*aq)

zeta_c4 = (2*kb*(1+alpha_i))/(1+2*kb*(1+alpha_i))

%slopes
% if ~ishandle(2)
%     f2 = figure('Position',[400,100,700,500]);
% else
%     figure(f2);
% end
% f2 = figure('Position',[400,100,700,500]);
% hold on
% % xlim([0,90])
% xlabel('$\phi$')
% ylabel('$\partial E/\partial \phi$')
% lines = ["-", ":", ":", "--", ":", "--"];
% colours = ['k', 'b', 'r', 'r', 'g', 'g'];
% phi_diffs = rad2deg((phi_vals(1:end-1)+phi_vals(2:end))/2);
% for ii=1:6
%     slopes(ii,:) = diff(E_all_linear(ii,:))./diff(rad2deg(phi_vals));
%     plot(phi_diffs,slopes(ii,:), ...
%         strcat(colours(ii),lines(ii)))
% end
% min_slope = min(slopes(1,:));
% ylim([-2*abs(min_slope), 2*abs(min_slope)])
% annotation('textbox', [0.5,0.7,0.4,0.2], 'String',...
%     [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
%     sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
%     sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
%     sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
%     sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta),...
%     sprintf('$\\kappa = %0.2g k_\\mathrm{D} R^2$ \n', kappa/(kD*R^2))], ...
%     'interpreter', 'latex','FontSize',16,...
%     'FitBoxToText','on','LineStyle','none', 'Color','b')
% legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
%     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
%     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
%     'location', 'nw')


% % areas
% figure();
% hold on
% xlim([0,90]);
% xlabel('$\phi$')
% yyaxis left
% ylabel('Area')
% plot(rad2deg(phi_vals), S_B_vals, '--', 'displayname', '$S_B$')
% plot(rad2deg(phi_vals), S_A_vals+S_B_vals, '-', 'displayname', '$S_\mathrm{total}$')
% yyaxis right
% ylabel('Area')
% plot(rad2deg(phi_vals), S_A_vals, '--', 'displayname', '$S_A$')
% legend
% 
% % stretches
% figure();
% hold on
% % xlim([0,90]);
% xlabel('$\phi$')
% ylabel('$\alpha$')
% plot(rad2deg(phi_vals), alpha_B_vals_linear, '-', 'displayname', '$\alpha_B$')
% plot(rad2deg(phi_vals), alpha_A_vals_linear, '-', 'displayname', '$\alpha_A$')
% plot(rad2deg(phi_vals), alpha_B_vals_toroid, '--', 'displayname', '$\alpha_B$')
% plot(rad2deg(phi_vals), alpha_A_vals_toroid, '--', 'displayname', '$\alpha_A$')
% plot(rad2deg(phi_vals), alpha_B_vals_nonlinear, ':', 'displayname', '$\alpha_B$')
% plot(rad2deg(phi_vals), alpha_A_vals_nonlinear, ':', 'displayname', '$\alpha_A$')
% legend

function [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_linear_minimum(constants, inputs)

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
    h_phi_init = inputs(3);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = 0;
    S_i = d^2;

    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;
    h_phi = h_phi_init;

    options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp');

    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init, h_phi_init],...
        [],[],[],[],[-0.1,-0.1,-Inf],[0.1, 0.1, Inf], ...
        @constraint, options);

    function f = objective(y_obj)

        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        h_phi = y_obj(3);
        
        Sigma = kD*alpha_B;
        r_phi = sin(phi)*R;
        r = linspace(r_phi, d/2,N);
        
        [~,~,S_B, ~,  ~, lap_h] = free_shape_linear_fixed_h(...
            r, r_phi, d, phi, kappa, Sigma, h_phi);

        % stretching, adhesion and bending energy
        f = (epsilon*n0*S_A./(1+alpha_A) ...
          + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
          + kappa/2*2*pi*trapz(r, r.*lap_h.^2) + 4*pi*kappa*(1-cos(phi)));

    end

    function [c,ceq]= constraint(~)
    
        c(1) = 0;
        ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
    
    end

end

% function f = free_linear(y, const)
%     % function of free energy to minimise
% 
%     epsilon = const(1);
%     n0 = const(2);
%     d = const(3);
%     R = const(4);
%     kD = const(5);
%     kappa = const(6);
%     alpha_i = const(7);
%     N = const(8);
%     phi = const(9);
%     
%     alpha_A = y(1);
%     alpha_B = y(2);
%     h_phi = y(3);
%     
%     Sigma = kD*alpha_B;
%     lambda = sqrt(kappa/Sigma);
%     r_phi = sin(phi)*R;
%     r = linspace(r_phi, d/2,N);
%     
%     [~,~,S_B, ~, lap_h, ~] = free_shape_linear_fixed_h(...
%         r, r_phi, d, phi, kappa, Sigma, h_phi);
%     
% %     hderiv = C(1)./r+C(3)/lambda*besseli(1,r/lambda)-C(4)/lambda*besselk(1,r/lambda);
% %     lap_h = C(3)/lambda^2*besseli(0,r/lambda)+C(4)/lambda^2*besselk(0,r/lambda);
%     S_A = 2*pi*R^2*(1-cos(phi));
%     
%     % stretching, adhesion and bending energy
%     f = epsilon*n0*S_A./(1+alpha_A) ...
%       + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
%       + kappa/2*2*pi*trapz(r, r.*lap_h.^2) + 4*pi*kappa*(1-cos(phi));
% 
% end
% 
% function [c,ceq] = area_con_linear(y, const)
% 
%     epsilon = const(1);
%     n0 = const(2);
%     d = const(3);
%     R = const(4);
%     kD = const(5);
%     kappa = const(6);
%     alpha_i = const(7);
%     N = const(8);
%     phi = const(9);
%     
%     alpha_A = y(1);
%     alpha_B = y(2);
%     h_phi = y(3);
%     
%     Sigma = kD*alpha_B;
%     r_phi = sin(phi)*R;
%     
%     r = linspace(r_phi, d/2,N);
%     [~,~,S_B, ~, ~, ~] = free_shape_linear_fixed_h(...
%         r, r_phi, d, phi, kappa, Sigma, h_phi);
%     
%     S_i = d^2;
%     S_A = 2*pi*R^2*(1-cos(phi));
% 
%     c(1) = 0;
% 
%     ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
% 
% end

function [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_nonlinear_minimum(constants, inputs)

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

    options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp');
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init],...
        [],[],[],[],[-0.1,eps(1)],[0.1, 0.1], ...
        @constraint, options);

    function f = objective(y_obj)

        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        
        Sigma = kD*alpha_B;
        lambda = sqrt(kappa/Sigma);

        psi_dot_init = -sqrt(2*(1-cos(phi)))/lambda;
        
        out_shape = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, psi_dot_init);
        
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = out_shape.ye(8)*2*pi + pi*((d/2)^2-out_shape.ye(2)^2) + d^2*(1-pi/4);
        
        % stretching, adhesion and bending energy
        f = epsilon*n0*S_A./(1+alpha_A) ...
          + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
          + out_shape.ye(7)*kappa*pi + 4*pi*kappa*(1-cos(phi));

    end

    function [c,ceq]= constraint(~)
    
        c(1) = 0;
        ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
    
    end

end

% function f = free_nonlinear(y, const)
%     % function of free energy to minimise
% 
%     epsilon = const(1);
%     n0 = const(2);
%     d = const(3);
%     R = const(4);
%     kD = const(5);
%     kappa = const(6);
%     alpha_i = const(7);
%     N = const(8);
%     phi = const(9);
%     psi_dot_init = const(10);
%     
%     alpha_A = y(1);
%     alpha_B = y(2);
%     
%     Sigma = kD*alpha_B;
%     
%     out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, psi_dot_init);
%     
%     S_A = 2*pi*R^2*(1-cos(phi));
%     S_B = out.ye(8)*2*pi + pi*((d/2)^2-out.ye(2)^2) + d^2*(1-pi/4);
%     
%     % stretching, adhesion and bending energy
%     f = epsilon*n0*S_A./(1+alpha_A) ...
%       + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
%       + out.ye(7)*kappa*pi + 4*pi*kappa*(1-cos(phi));
% 
%     if out.ye(1)>1e-3
%         f = f*1e5;
%     end
% 
% end
% 
% function [c,ceq] = area_con_nonlinear(y, const)
% 
%     epsilon = const(1);
%     n0 = const(2);
%     d = const(3);
%     R = const(4);
%     kD = const(5);
%     kappa = const(6);
%     alpha_i = const(7);
%     N = const(8);
%     phi = const(9);
%     psi_dot_init = const(10);
%     
%     alpha_A = y(1);
%     alpha_B = y(2);
%     
%     Sigma = kD*alpha_B;
%     
%     out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, psi_dot_init);
%     
%     S_i = d^2;
%     S_A = 2*pi*R^2*(1-cos(phi));
%     S_B = out.ye(8)*2*pi + pi*((d/2)^2-out.ye(2)^2) + d^2*(1-pi/4);
% 
%     c(1) = 0;
% 
%     ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
% 
% end

function [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_toroid_minimum(constants, inputs)

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
    rho_init = inputs(3);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = 0;
    S_i = d^2;

    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;
    rho = rho_init;

    options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp');
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init, rho_init],...
        [],[],[],[],[-0.1,-0.1, 0],[0.1, 0.1, 5*d], ...
        @constraint, options);

    function f = objective(y_obj)
        
        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        rho = y_obj(3); %toroid radius
    
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

end

% function f = free_toroid(y, const)
%     % function of free energy to minimise
% 
%     epsilon = const(1);
%     alpha_i = const(2);
%     d = const(3);
%     R = const(4);
%     phi = const(5);
%     kappa = const(6);
%     kD = const(7);
% 
%     alpha_A = y(1);
%     alpha_B = y(2);
%     rho = y(3); %toroid radius
% 
%     delta = sin(phi)*(R+rho);
% 
%     S_A = 2*pi*R^2*(1-cos(phi));
%     S_B = d^2 - pi*delta^2 + 2*pi*rho*(delta*phi-rho*(1-cos(phi)));
% 
%     E_adhesion = epsilon*S_A./(1+alpha_A);
%     E_stretch_A = kD/2*alpha_A.^2*S_A./(1+alpha_A);
%     E_stretch_B = kD/2*alpha_B.^2*S_B./(1+alpha_B);
%     E_bend_A = 4*pi*kappa*(1-cos(phi));
%     a = delta/rho;
%     if a>1
%         E_bend_B = 2*pi*kappa*(a^2/sqrt(a^2-1)*...
%             (acot(sqrt(a^2-1))+atan((a*tan(phi/2)-1)/sqrt(a^2-1)))-2*(1-cos(phi)));
%     else
%         E_bend_B = 2*pi*kappa*(a^2/sqrt(1-a^2)*...
%             (-acoth(sqrt(1-a^2))+atanh((1-a*tan(phi/2))/sqrt(1-a^2)))-2*(1-cos(phi)));
%     end
% 
%     f = E_adhesion+E_stretch_A+E_stretch_B+E_bend_A+E_bend_B;
% 
% end
% 
% function [c,ceq] = area_con_toroid(y, const)
% 
%     epsilon = const(1);
%     alpha_i = const(2);
%     d = const(3);
%     R = const(4);
%     phi = const(5);
%     kappa = const(6);
%     kD = const(7);
% 
%     alpha_A = y(1);
%     alpha_B = y(2);
%     rho = y(3);
% 
%     delta = sin(phi)*(R+rho);
% 
%     S_i = d^2;
%     S_A = 2*pi*R^2*(1-cos(phi));
%     S_B = d^2 - pi*delta^2 + 2*pi*rho*(delta*phi-rho*(1-cos(phi)));
% 
%     c(1) = 0;
% 
%     ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
% 
% end
