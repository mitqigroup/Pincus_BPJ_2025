function [r,h,out] = ...
    free_shape_nonlinear_fixed_h(R, d, phi, kappa, Sigma, h_phi,psi_dot_init, p_r_init)
% now solving the nonlinear shape equation using shooting method
    
    lambda = sqrt(kappa/Sigma);
    const(1) = lambda;
    const(2) = d;
    
    r_phi = R*sin(phi);
    psi_init = phi;
%     p_psi_init = 2*r_phi*(psi_dot_init+sin(psi_init)/r_phi);
%     z = 1-cos(phi);
%     p_r_init = sqrt(z*(2-z))/(1-z)*(1/R+2*z*R/lambda^2-R*psi_dot_init^2);

    min_constants = [lambda, R, phi, d, h_phi];

    options = optimoptions('fminunc', 'MaxFunctionEvaluations',1000)

%     outputs = fminsearch(@(inp) get_shape(inp, min_constants), [psi_dot_init, p_r_init]);
    outputs = fminunc(@(inp) get_shape(inp, min_constants), [psi_dot_init, p_r_init],...
        options);

    psi_dot_init_final = outputs(1);
    p_r_init_final = outputs(2);
    p_psi_init_final = 2*r_phi*(psi_dot_init_final+sin(psi_init)/r_phi);

    init_vals = [psi_init, r_phi, R*sin(phi), p_psi_init_final, p_r_init_final, 0,0];
    event_func = @(s,y) myEventFcn(s,y,const);
%     options = odeset('Events', @myEventFcn);
    options = odeset('Events', event_func);
    out = ode45(@(s, y) hamilton(s, y,const),linspace(0,10,1000),init_vals,...
        options);

    solution = deval(out, linspace(0,out.xe, 1000));
    r = solution(2,:);
    h = solution(3,:)-R*sin(phi)+h_phi;
%     r = out.y(2,:);
%     h = out.y(3,:)-R*sin(phi)+h_phi;

%     init_vals = [psi_init, r_phi, h_phi, p_psi_init, p_r_init, 0]
%     out = ode45(@(s, y) hamilton(s, y,const),[0,10],init_vals);
%     r = out.y(2,:);
%     h = out.y(3,:);

end

function derivs = hamilton(s, y, const)
    psi = y(1);
    r = y(2);
    h = y(3);
    p_psi = y(4);
    p_r = y(5);
    p_h = y(6);
    E_bend = y(7);

    lambda = const(1);
    d = const(2);

    f(1) = p_psi/(2*r)-sin(psi)/r;
    f(2) = cos(psi);
    f(3) = sin(psi);
    f(4) = cos(psi)*(p_psi/r-p_h)...
        +sin(psi)*(2*r/lambda^2+p_r);
    f(5) = p_psi/r*(p_psi/(4*r)-sin(psi)/r)...
        +2/lambda^2*(1-cos(psi));
    f(6) = 0;
    f(7) = p_psi^2/(2*r);

    derivs = f';
end

function [value,isterminal,direction] = myEventFcn(s,y,const)
    psi = y(1);
    r = y(2);
    h = y(3);
    p_psi = y(4);
    p_r = y(5);
    p_h = y(6);

    lambda = const(1);
    d = const(2);

    value = y(2)-d/2;
    isterminal = 1;
    direction = 0;
end

function diffs = get_shape(inputs, constants)

    psi_dot_init = inputs(1);
    p_r_init = inputs(2);

    lambda = constants(1);
    R = constants(2);
    phi = constants(3);
    d = constants(4);
    h_phi = constants(5);

    const(1) = lambda;
    const(2) = d;
    
    r_phi = R*sin(phi);
    psi_init = phi;
    p_psi_init = 2*r_phi*(psi_dot_init+sin(psi_init)/r_phi);

    init_vals = [psi_init, r_phi, R*sin(phi), p_psi_init, p_r_init, 0,0];
    event_func = @(s,y) myEventFcn(s,y,const);
%     options = odeset('Events', @myEventFcn);
    options = odeset('Events', event_func);
    out = ode45(@(s, y) hamilton(s, y,const),[0,100],init_vals,...
        options);
    h_end = out.y(3,end)-R*sin(phi)+h_phi;
    psi_end = out.y(1,end);

    diffs = (h_end^2+psi_end^2);
end