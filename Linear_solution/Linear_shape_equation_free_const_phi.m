% clear variables

% constants
R = 0.3;
sigma = 3.5e-3;
d = sqrt(R^2/sigma);    % um
% d = 20;    % um
phi = pi/12;
kD  = 300/10^12*1e9;    % picoJ/um^2
zeta = 7.43e-6;            % dimensionless
epsilon = -zeta*kD;     % picoJ/um^2
n0 = 1;                 % fraction
kappa = 1e-19*1e12;     % picoJ
kappa_bar = kappa/(kD*R^2);
% gamma = 0.5*R;            % um
% epsilon = -kappa/gamma^2;
% kappa = 1e-17*1e12;     % picoJ
% zeta = -epsilon*n0/kD;   % dimensionless
% alpha_i = sqrt(-zeta/2);        % fraction
alpha_i = 0.003;
% alpha_i = 0;
% ten_param=-d*epsilon*n0/sqrt(kD*kappa)
N = 3e3;                % number of points in quadrature

savedata = 0;
plotfigs = 0;
filename = '300nm_ai1kap19zeta002.mat';

phi_vals = flip(deg2rad(linspace(0.001,0.1,50)));

% phi_vals = deg2rad(28);

E_all = zeros(6,length(phi_vals));

tic
for ii = 1:length(phi_vals)
phi = phi_vals(ii);
rad2deg(phi)

if ii==1
    alpha_A_init = (alpha_i-zeta)/(1+zeta);
    alpha_B_init = alpha_i*2;
else
    alpha_A_init = alpha_A;
    alpha_B_init = alpha_B;
end

%% use initial stretch to solve for minimum attachment
% options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e3);
% options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp');
% const = [epsilon, n0, d, R, kD, kappa, alpha_i, N, phi];
% [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
%     fmincon(@(y) stretch_bend_min(y, const),[alpha_A_init, alpha_B_init, h_phi_init],...
%     [],[],[],[],[-0.1,-0.1,-Inf],[0.1, 0.1, Inf], ...
%     @(y) lipid_con_bend(y,const), options);

const = [epsilon, n0, d, R, kD, kappa, alpha_i, phi];
[out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_linear_minimum(const, [alpha_A_init, alpha_B_init]);

alpha_A = out(1);
alpha_B = out(2);

lambda_stretch = lam_vals.eqnonlin;

% get the shape of the free region
Sigma = kD*alpha_B;
lambda = sqrt(kappa/Sigma);
r_phi = sin(phi)*R;
r = linspace(r_phi, d/2,N);

[C, delA, Ebend] = free_shape_linear_free_h(R, phi, kappa, Sigma);

S_A = 2*pi*R^2*(1-cos(phi));
S_B = delA+(d^2-pi*r_phi^2);

E_adhesion = epsilon*n0*S_A./(1+alpha_A);
E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
E_bend_B = Ebend;
E_bend_A = 4*pi*kappa*(1-cos(phi));

h = C(1)+C(2)*besselk(0,r/lambda);

E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;

initial_E = alpha_i^2*kD*d^2/(1+alpha_i)/2;

E_all(1,ii) = E-initial_E;
E_all(2,ii) = E_adhesion;
E_all(3,ii) = E_stretch_A;
E_all(4,ii) = E_stretch_B-initial_E;
E_all(5,ii) = E_bend_A;
E_all(6,ii) = E_bend_B;

alpha_A_vals(ii) = alpha_A;
alpha_B_vals(ii) = alpha_B;
S_A_vals(ii) = S_A;
S_B_vals(ii) = S_B;
Sigma_vals(ii) = Sigma;

%% plot of shape and nanoparticle
% subplot(1,2,1);
if plotfigs
figure('Position',[400,100,800,600]);
hold on
% axis equal
xlabel('$r$')
ylabel('$h$')
plot(r, h, 'displayname', 'free surface');
t = linspace(-pi/2,pi/2,1000);
x = cos(t)*R;
% y = sin(t)*R+(R*cos(phi)+h(1));
y = sin(t)*R+R*cos(phi);
plot(x,y, 'displayname', 'microbead')

end

end
toc

if savedata
    save(filename);
end

% if ~ishandle(1)
%     f1 = figure('Position',[400,100,700,500]);
% else
%     figure(f1);
% end
% f1 = figure('Position',[400,100,700,500]);
hold on
% xlim([0,90])
xlabel('$\phi$')
ylabel('$\Delta E$')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
for ii=[1]
    plot(rad2deg(phi_vals), E_all(ii,:),'-',...
        'displayname', sprintf('$\\zeta = %d$', zeta))
end

zeta_c = (2*kappa_bar*(1+alpha_i))/(1-alpha_i-1/2*alpha_i^2+2*kappa_bar*(1+alpha_i))

kb = kappa_bar;
% ai = alpha_i;
% 
% aq = (20-12*ai-4*ai^2+24*kb*(1+ai));
% bq = -(-8+8*ai+4*ai^2-16*kb*(1+ai));
% cq = (16*kb*(1+ai));
% 
% zeta_c2 = (-bq-sqrt(bq^2-4*aq*cq))/(2*aq)
% zeta_c3 = (-bq+sqrt(bq^2-4*aq*cq))/(2*aq)
% 
% aq = (20+24*kb*(1+ai));
% bq = -(-8-16*kb*(1+ai));
% cq = (16*kb*(1+ai));
% zeta_c5 = (-bq+sqrt(bq^2-4*aq*cq))/(2*aq)
% 
% zeta_c4 = (2*kb*(1+alpha_i))/(1+2*kb*(1+alpha_i))

function [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
    get_linear_minimum(constants, inputs)

    epsilon = constants(1);
    n0 = constants(2);
    d = constants(3);
    R = constants(4);
    kD = constants(5);
    kappa = constants(6);
    alpha_i = constants(7);
    phi = constants(8);

    grad = [];
    hessian = [];
    
    alpha_A_init = inputs(1);
    alpha_B_init = inputs(2);

    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;
    
    Sigma = kD*alpha_B;
    r_phi = sin(phi)*R;
    lambda = sqrt(kappa/Sigma);
    lamt = lambda*sqrt(alpha_B);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = (d^2-pi*r_phi^2);
    S_i = d^2;

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%         'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
%         'SpecifyObjectiveGradient', true, 'CheckGradients', true);

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%         'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
%         'SpecifyObjectiveGradient', true);

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%         'OptimalityTolerance', 1e-14, 'ConstraintTolerance', 1e-14, 'StepTolerance', 1e-14,...
%         'SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient', true,...
%         'Diagnostics', 'on', 'Display', 'iter-detailed');

    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
        'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
        'SpecifyObjectiveGradient', false);

    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@objective,[alpha_A_init, alpha_B_init],...
        [],[],[],[],[-0.1,-0.1],[0.1, 0.1], ...
        @constraint, options);
% 
%     clear options
%     options.verify_level = 1;
%     options.printfile = 'snsolve_print.txt';
%     options.scale_option = 2;
%     options.major_feasibility_tolerance = 1e-12;
%     options.major_optimality_tolerance = 1e-12;
% %     options.function_precision = 1e-8;
% %     options.major_optimality_tolerance = 1e-4;
%     [out,fval,exitflag,output,lam_vals,states] = ...
%         snsolve(@objective,[alpha_A_init, alpha_B_init], [],[],[],[],...
%         [-0.1,-0.1]',[0.1, 0.1]', ...
%         @constraint, options);

    function [f,g] = objective(y_obj)

        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        
        Sigma = kD*alpha_B;
        r_phi = sin(phi)*R;
        lambda = sqrt(kappa/Sigma);
        lamt = lambda*sqrt(alpha_B);

        [~, delA, Ebend] =...
            free_shape_linear_free_h(R, phi, kappa, Sigma);
    
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = delA+(d^2-pi*r_phi^2);

        alpha_A_gradient = (kD/2*S_A*alpha_A*(alpha_A+2)-epsilon*n0*S_A)...
            /(alpha_A+1)^2;

%         alpha_A_gradient = -epsilon*n0*S_A./(1+alpha_A)^2+kD*alpha_A*S_A/(1+alpha_A)...
%             -kD*1/2*alpha_A^2*S_A/(1+alpha_A)^2

        parts(1) = -kD*1/2*d^2*alpha_B^2/(1+alpha_B)^2;
        parts(2) = kD*d^2*alpha_B/(1+alpha_B);
        parts(3) = kD/2*pi*r_phi^2*alpha_B^2/(1+alpha_B)^2;
        parts(4) = -kD*pi*r_phi^2*alpha_B/(1+alpha_B);
        parts(5) = -kD*pi*r_phi^2*alpha_B^2*tan(phi)^2/(1+alpha_B)^2/4;
        parts(6) = +kD*3/4*r_phi^2*alpha_B*tan(phi)^2/(1+alpha_B)^2;
        parts(7) = -pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2);
        parts(8) = +(pi*besselk(0,r_phi/lambda,1)* ...
        (-4*kD*r_phi*sqrt(alpha_B)*lamt^4 ...
         -2*kD*r_phi*alpha_B^(3/2)*lamt^4 ...
         -2*r_phi^3*sqrt(alpha_B)*kappa ...
         -4*r_phi^3*alpha_B^(3/2)*kappa ...
         -2*r_phi^3*alpha_B^(5/2)*kappa ...
         +kD*r_phi^3 *alpha_B^(3/2)*lamt^2 ...
         +kD*r_phi^3*alpha_B^(5/2)*lamt^2)*tan(phi)^2)...
         /(4*(1+alpha_B)^2*lamt^3*besselk(1,r_phi/lambda,1));
        parts(9) = (besselk(0,r_phi/lambda,1)^2)* ...
         ((kD*pi*r_phi^2*tan(phi)^2*alpha_B^2)/(4*(1+alpha_B)^2) ...
         -(kD*5*pi*r_phi^2*tan(phi)^2*alpha_B)/(8*(1+alpha_B)) ...
         +pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2) ...
         +(pi*r_phi^2*(4*kappa+4*kD*alpha_B*kappa-3*kD*alpha_B*lamt^2)*tan(phi)^2 ...
         /(8*(1+alpha_B)*lamt^2))...
         )/(besselk(1,r_phi/lambda,1))^2;
        parts(10) = (besselk(0,r_phi/lambda,1)^3* ...
         (pi*r_phi^3*sqrt(alpha_B)*kappa*tan(phi)^2/(4*lamt^3)  ...
         -kD*pi*r_phi^3*alpha_B^(3/2)*tan(phi)^2/(8*(1+alpha_B)*lamt) ...
         +pi*r_phi^3*tan(phi)^2*(2*sqrt(alpha_B)*kappa+2*alpha_B^(3/2) ...
         *kappa-kD*alpha_B^(3/2)*lamt^2)/(8*(1+alpha_B)*lamt^3))) ...
         /(besselk(1,r_phi/lambda,1)^3);

        alpha_B_gradient = (...
        -kD*1/2*d^2*alpha_B^2/(1+alpha_B)^2 ...
        +kD*d^2*alpha_B/(1+alpha_B) ...
        +kD/2*pi*r_phi^2*alpha_B^2/(1+alpha_B)^2 ... 
        -kD*pi*r_phi^2*alpha_B/(1+alpha_B) ...
        -kD*pi*r_phi^2*alpha_B^2*tan(phi)^2/(1+alpha_B)^2/4 ...
        +kD*3/4*r_phi^2*alpha_B*tan(phi)^2/(1+alpha_B)^2 ...
        -pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2)...
        +(pi*besselk(0,r_phi/lambda,1)* ...
        (-4*kD*r_phi*sqrt(alpha_B)*lamt^4 ...
         -2*kD*r_phi*alpha_B^(3/2)*lamt^4 ...
         -2*r_phi^3*sqrt(alpha_B)*kappa ...
         -4*r_phi^3*alpha_B^(3/2)*kappa ...
         -2*r_phi^3*alpha_B^(5/2)*kappa ...
         +kD*r_phi^3 *alpha_B^(3/2)*lamt^2 ...
         +kD*r_phi^3*alpha_B^(5/2)*lamt^2)*tan(phi)^2)...
         /(4*(1+alpha_B)^2*lamt^3*besselk(1,r_phi/lambda,1)) ...
        + (besselk(0,r_phi/lambda,1)^2)* ...
         ((kD*pi*r_phi^2*tan(phi)^2*alpha_B^2)/(4*(1+alpha_B)^2) ...
         -(kD*5*pi*r_phi^2*tan(phi)^2*alpha_B)/(8*(1+alpha_B)) ...
         +pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2) ...
         +(pi*r_phi^2*(4*kappa+4*kD*alpha_B*kappa-3*kD*alpha_B*lamt^2)*tan(phi)^2 ...
         /(8*(1+alpha_B)*lamt^2))...
         )/(besselk(1,r_phi/lambda,1))^2 ...
        + (besselk(0,r_phi/lambda,1)^3* ...
         (pi*r_phi^3*sqrt(alpha_B)*kappa*tan(phi)^2/(4*lamt^3)  ...
         -kD*pi*r_phi^3*alpha_B^(3/2)*tan(phi)^2/(8*(1+alpha_B)*lamt) ...
         +pi*r_phi^3*tan(phi)^2*(2*sqrt(alpha_B)*kappa+2*alpha_B^(3/2) ...
         *kappa-kD*alpha_B^(3/2)*lamt^2)/(8*(1+alpha_B)*lamt^3))) ...
         /(besselk(1,r_phi/lambda,1)^3));

        g = [alpha_A_gradient, alpha_B_gradient]';

        E_adhesion = epsilon*n0*S_A./(1+alpha_A);
        E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
        E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
        E_bend_B = Ebend;
        E_bend_A = 4*pi*kappa*(1-cos(phi));

        % stretching, adhesion and bending energy
        f = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;

    end

    function [c,ceq, gc, gceq]= constraint(~)
    
        c = [];
        ceq = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

        gc = [];

%         alpha_A
%         alpha_B
%         S_A
%         S_B
%         ceq

%         parts(1) = 

        alpha_A_gradient = -S_A/(1+alpha_A)^2;
%         alpha_B_gradient_old = 1/(1+alpha_B)^2*( ...
%             -d^2+pi*r_phi^2-pi*r_phi^2*tan(phi)^2/2 ...
%             +(pi*besselk(0,r_phi/lambda,1))/(2*sqrt(alpha_B)*lamt*besselk(1,r_phi/lambda,1)) ...
%             *(2*r_phi*lamt^2+r_phi^3+r_phi^3*alpha_B)*tan(phi)^2) ...
%             + 1/(1+alpha_B)*(pi*r_phi^2*tan(phi)^2/2 ...
%             -pi*r_phi^3*tan(phi)^2*besselk(0,r_phi/lambda,1)^3 ...
%             /besselk(1,r_phi/lambda,1)^3/(2*sqrt(alpha_B))) ...
%             +besselk(0,r_phi/lambda,1)^2/besselk(1,r_phi/lambda,1)^2 ...
%             *(pi*r_phi^2*tan(phi)/(2*(1+alpha_B)^2)-pi*r_phi^2*tan(phi)/(alpha_B*(1+alpha_B)));

        alpha_B_gradient = (-1).*d.^2.*(1+alpha_B).^(-2)+pi.*R.^2.*(1+alpha_B).^(-2).*sin(phi).^2+(-1/2) ...
          .*pi.*R.^2.*(1+alpha_B).^(-2).*sin(phi).^2.*tan(phi).^2+(1/2).*pi.*R.^2.* ...
          alpha_B.^(-1).*(1+alpha_B).^(-1).*sin(phi).^2.*tan(phi).^2+(-1/2).*pi.*R.^3.* ...
          alpha_B.^(-1/2).*(1+alpha_B).^(-1).*lamt.^(-1).*besselk(0,R.*alpha_B.^(1/2).*lamt.^( ...
          -1).*sin(phi),1).^3.*besselk(1,R.*alpha_B.^(1/2).*lamt.^(-1).*sin(phi),1).^(-3).* ...
          sin(phi).^3.*tan(phi).^2+(1/2).*pi.*alpha_B.^(-1/2).*(1+alpha_B).^(-2).*lamt.^(-1) ...
          .*besselk(0,R.*alpha_B.^(1/2).*lamt.^(-1).*sin(phi),1).*besselk(1,R.*alpha_B.^( ...
          1/2).*lamt.^(-1).*sin(phi),1).^(-1).*(2.*R.*lamt.^2.*sin(phi)+R.^3.*sin(phi) ...
          .^3+R.^3.*alpha_B.*sin(phi).^3).*tan(phi).^2+besselk(0,R.*alpha_B.^(1/2).*lamt.^( ...
          -1).*sin(phi),1).^2.*besselk(1,R.*alpha_B.^(1/2).*lamt.^(-1).*sin(phi),1).^(-2).* ...
          ((1/2).*pi.*R.^2.*(1+alpha_B).^(-2).*sin(phi).^2.*tan(phi).^2+(-1).*pi.* ...
          R.^2.*alpha_B.^(-1).*(1+alpha_B).^(-1).*sin(phi).^2.*tan(phi).^2);


        gceq = [alpha_A_gradient,  alpha_B_gradient]';


    end

end