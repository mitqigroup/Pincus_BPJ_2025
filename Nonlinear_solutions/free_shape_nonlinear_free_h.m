function out = ...
    free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, debug_on)
% solve nonlinear shape equation using hunt and bracket

    if Sigma<0
        Sigma = -Sigma;
    end

    warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
    lambda = sqrt(kappa/Sigma);
    sigma = R^2/lambda^2;
    constants = [lambda, R, phi, d];
    z = 1-cos(phi);
    terminate_condition = 1e-10;

%     debug_on = false;
%     debug_on = true;

    if debug_on
        % psi_0_dot = linspace(-sqrt(sigma)/R-1,1/R,200);
        psi_0_dot = linspace(-4*sqrt(sigma)/R-2/R,2/R,200);
        % psi_0_dot = linspace(-1.16,-1.13,500);
        % psi_0_dot = linspace(-1e6,1,500);
        s_end = zeros(size(psi_0_dot));
        correct = zeros(size(psi_0_dot));
        correct2 = zeros(size(psi_0_dot));
        for ii = 1:length(psi_0_dot)
            [~,out] = get_shape(psi_0_dot(ii), constants, 1);
            s_end(ii) = out.x(end);
            correct(ii) = isempty(out.xe)||out.ye(1)>0.01;
            correct2(ii) = isempty(out.xe);
        end
    
        figure('Position',[400,100,800,600]);
        hold on
        plot(psi_0_dot, s_end, '-')
        xlabel('$\dot{psi_0}$')
        ylabel('$s_\mathrm{max}$')
        yyaxis right
        ylim([-0.05,1.05])
        plot(psi_0_dot, correct)
        ylabel('$\phi =\pi$ True/False')
        annotation("textbox","String",{sprintf('$\\phi = %0.3g^{\\circ}$', rad2deg(phi)),...
            sprintf('$R = %0.2g$', R), sprintf('$\\bar{\\sigma} = %0.2g$', sigma)});
        % xlim([-20*sqrt(sigma)/R-1,1/R])
        yyaxis left
    end
    
    psi_0_dot_upper = 1/R;
    psi_0_dot_lower = -4*sqrt(sigma)/R-2/R;
    
    N = 1000;
    psi_0_dot_iter = zeros(1,N);
    s_end_iter = zeros(1,N);
    physical_solution = false(1,N);
    can_find_solution = true;
    
    psi_0_dot_iter_up = zeros(1,N);
    s_end_iter_up = zeros(1,N);
    physical_solution_up = false(1,N);
    
    consider_up_solution = false;
    greater_than_d_solution = false;
    consider_up_solution_greater_than_d = false;
    consider_down_solution_greater_than_d = false;
    one_on_R_valid = false;
    
    options = optimoptions('fmincon', 'algorithm','sqp', 'MaxFunctionEvaluations',1000,...
        'OptimalityTolerance', 1e-9, 'ConstraintTolerance', 1e-9, 'StepTolerance', 1e-9, ...
        'Display','none');
    
    if z<1
        psi_0_dot_iter(1) = psi_0_dot_lower;
        alpha = sqrt(sigma)/(10*R);
        counter = 1;
        crossed = false;
        while true
            % start by hunting upwards, until we cross
            if ~crossed
                [~,out] = get_shape(psi_0_dot_iter(counter), constants, 1);
                if isempty(out.xe)||out.ye(1)>0.01
                    % we've crossed, bracket
                    crossed = true;
                    psi_0_dot_upper = psi_0_dot_iter(counter);
                    psi_0_dot_lower = psi_0_dot_iter(counter-1);
                else
                    alpha = alpha*2;
                    psi_0_dot_iter(counter+1) = psi_0_dot_iter(counter) + alpha; 
                    if psi_0_dot_iter(counter+1)>1/R
                        % we've hit the unphysical limit, bracket
                        crossed = true;
                        psi_0_dot_upper = 1/R;
                        psi_0_dot_lower = psi_0_dot_iter(counter);
                    end
                end
            else
                %if we have crossed, now bisect
                psi_0_dot_iter(counter) = (psi_0_dot_upper+psi_0_dot_lower)/2;
                [~,out] = get_shape(psi_0_dot_iter(counter), constants, 1);
                if isempty(out.xe)||out.ye(1)>0.01
                    % we're above the limit, current psi_0_dot_dot is new upper
                    psi_0_dot_upper = psi_0_dot_iter(counter);
                else
                    % we're below the limit, current psi_0_dot_dot is new lower
                    psi_0_dot_lower = psi_0_dot_iter(counter);
                end
                alpha = abs(psi_0_dot_upper-psi_0_dot_lower);
                % separate check - if we have a physical point with r_end>d/2,
                % we instead solve directly using fmincon
                if ~(isempty(out.xe)||out.ye(1)>0.01)&&(out.ye(2)>d/2)
                    outputs_up = fmincon(@(inp) get_shape(inp, constants, 2), psi_0_dot_iter(counter),...
                                [],[],[],[],-4*sqrt(sigma)/R-2/R,psi_0_dot_iter(counter),[],options);
                    greater_than_d_solution = true;
                    consider_up_solution_greater_than_d = true;
                    break;
                end
            end
    
            s_end_iter(counter) = out.x(end);
            physical_solution(counter) = ~(isempty(out.xe)||out.ye(1)>0.01);
    
            if debug_on
                % for debugging
%                 fprintf('s iter = %0.2g\n',s_end_iter(counter))
%                 fprintf('psi 0 iter = %0.2g\n',psi_0_dot_iter(counter))
%                 fprintf('alpha = %0.2g\n',alpha)
                text(psi_0_dot_iter(counter), s_end_iter(counter), sprintf('%0.2g', counter))
                ylim([0,max(s_end_iter)])
            end
    
            if alpha<terminate_condition
                break
            else
                counter = counter + 1;
            end
        end
    else
    %     psi_0_dot_iter(1) = 1/(R);
    %     psi_0_dot_iter(1) = psi_0_dot(find(correct==0,1,'last'));
        psi_0_dot_iter(1) = -3/(10*R);
        % check that this first point is valid
        [~,out] = get_shape(psi_0_dot_iter(1), constants, 1);
        times_searched = 0;
        counter = 1;
        crossed = false;
        alpha = sqrt(sigma)/(10*R);
        while (isempty(out.xe)||out.ye(1)>0.01)
            % if not valid, start at 1/R and search down
            psi_0_dot_iter(1) = 1/R - alpha*times_searched/2;
            [~,out] = get_shape(psi_0_dot_iter(1), constants, 1);
            times_searched = times_searched + 1;
            if psi_0_dot_iter(1)<(-4*sqrt(sigma)/R-2/R)
                % at this point, just give up. I've got no idea what to do
                % here, maybe we can pick a close point in the phi-space? I
                % think for now just calculate at psi_0_dot = 1/R
                can_find_solution = false;
                psi_0_dot_iter(1) = 1/R;
                [~,out] = get_shape(psi_0_dot_iter(1), constants, 1);
                s_end_iter(counter) = out.x(end);
                physical_solution(counter) = ~(isempty(out.xe)||out.ye(1)>0.01);
                break
            end
    
        end
        while can_find_solution
            % start by hunting downwards, until we cross
            if ~crossed
                if counter>1
                    [~,out] = get_shape(psi_0_dot_iter(counter), constants, 1);
                end
                if isempty(out.xe)||out.ye(1)>0.01
                    % we've crossed, bracket
                    crossed = true;
                    psi_0_dot_lower = psi_0_dot_iter(counter);
                    psi_0_dot_upper = psi_0_dot_iter(counter-1);
                else
                    alpha = alpha*2;
                    psi_0_dot_iter(counter+1) = psi_0_dot_iter(counter) - alpha; 
                end
            else
                %if we have crossed, now bisect
                psi_0_dot_iter(counter) = (psi_0_dot_upper+psi_0_dot_lower)/2;
                [~,out] = get_shape(psi_0_dot_iter(counter), constants, 1);
                if isempty(out.xe)||out.ye(1)>0.01
                    % we're below the limit, current psi_0_dot is new upper
                    psi_0_dot_lower = psi_0_dot_iter(counter);
                else
                    % we're above the limit, current psi_0_dot is new lower
                    psi_0_dot_upper = psi_0_dot_iter(counter);
                end
                alpha = abs(psi_0_dot_upper-psi_0_dot_lower);
                % separate check - if we have a physical point with r_end>d/2,
                % we instead solve directly using fmincon
                if ~(isempty(out.xe)||out.ye(1)>0.01)&&(out.ye(2)>d/2)
                    outputs_down = fmincon(@(inp) get_shape(inp, constants, 2), psi_0_dot_iter(counter),...
                                [],[],[],[],psi_0_dot_iter(counter),psi_0_dot_iter(1),[],options);
                    greater_than_d_solution = true;
                    consider_down_solution_greater_than_d = true;
                    break;
                end
            end
    
            s_end_iter(counter) = out.x(end);
            physical_solution(counter) = ~(isempty(out.xe)||out.ye(1)>0.01);
    
            if debug_on
            % for debugging
%             fprintf('s iter = %0.2g\n',s_end_iter(counter))
%             fprintf('psi 0 iter = %0.2g\n',psi_0_dot_iter(counter))
%             fprintf('alpha = %0.2g\n',alpha)
            text(psi_0_dot_iter(counter), s_end_iter(counter), sprintf('%0.2g', counter))
            ylim([0,max(s_end_iter)])
            end
    
            if alpha<terminate_condition
                break
            else
                counter = counter + 1;
            end
        end
        % we also want to hunt up, just to be sure
        counter = 1;
        crossed = false;
        psi_0_dot_iter_up(1) = psi_0_dot_iter(1);
        alpha = sqrt(sigma)/(10*R);
        while can_find_solution
            % start by hunting upwards, until we cross or get past 1/R
            if ~crossed
                [~,out] = get_shape(psi_0_dot_iter_up(counter), constants, 1);
                if psi_0_dot_iter_up(counter)>1/R
                    % check that the 1/R solution is feasible, if it is, we
                    % should stop the search, if it isn't, we should bisect!
                    [~,out_test] = get_shape(1/R, constants, 1);
                    if isempty(out_test.xe)||out_test.ye(1)>0.01
                        % move to bisection
                        psi_0_dot_upper = 1/R;
                        psi_0_dot_lower = psi_0_dot_iter_up(counter-1);
                        crossed = true;
                    else
                        physical_solution_up(counter) = 0;
                        s_end_iter_up(counter) = 0;
                        break;
                    end
                end
                if isempty(out.xe)||out.ye(1)>0.01
                    % we've crossed, bracket
                    crossed = true;
                    psi_0_dot_upper = psi_0_dot_iter_up(counter);
                    psi_0_dot_lower = psi_0_dot_iter_up(counter-1);
                else
                    alpha = alpha*2;
                    psi_0_dot_iter_up(counter+1) = psi_0_dot_iter_up(counter) + alpha; 
                end
            else
                %if we have crossed, now bisect
                psi_0_dot_iter_up(counter) = (psi_0_dot_upper+psi_0_dot_lower)/2;
                [~,out] = get_shape(psi_0_dot_iter_up(counter), constants, 1);
                if isempty(out.xe)||out.ye(1)>0.01||(out.ye(1)>0.01&&out.ye(2)>d/2)
                    % we're below the limit, current psi_0_dot is new upper
                    psi_0_dot_upper = psi_0_dot_iter_up(counter);
                else
                    % we're above the limit, current psi_0_dot is new lower
                    psi_0_dot_lower = psi_0_dot_iter_up(counter);
                end
                alpha = abs(psi_0_dot_upper-psi_0_dot_lower);
                % separate check - if we have a physical point with r_end>d/2,
                % we instead solve directly using fmincon
                if ~(isempty(out.xe)||out.ye(1)>0.01)&&(out.ye(2)>d/2)
                    outputs_up = fmincon(@(inp) get_shape(inp, constants, 2), psi_0_dot_iter(counter),...
                                [],[],[],[],psi_0_dot_iter(1),psi_0_dot_iter(counter),[],options);
                    greater_than_d_solution = true;
                    consider_up_solution_greater_than_d = true;
                    break;
                end
            end
    
            s_end_iter_up(counter) = out.x(end);
            physical_solution_up(counter) = ~(isempty(out.xe)||out.ye(1)>0.01);
    
            if debug_on
                % for debugging
%                 fprintf('s iter = %0.2g\n',s_end_iter_up(counter))
%                 fprintf('psi 0 iter = %0.2g\n',psi_0_dot_iter_up(counter))
%                 fprintf('alpha = %0.2g\n',alpha)
                text(psi_0_dot_iter_up(counter), s_end_iter_up(counter), sprintf('%0.2g', counter), 'color', 'red')
    %             ylim([0,max(s_end_iter)])
            end
    
            if alpha<terminate_condition
                consider_up_solution = true;
                break
            else
                counter = counter + 1;
            end
        end
        % finally, we always want to check the 1/R point if we haven't already,
        % then hunt downwards
        if psi_0_dot_iter(1) ~= 1/R
            psi_0_dot_iter_1_R(1) = 1/R;
            % check that this first point is valid
            [~,out] = get_shape(psi_0_dot_iter_1_R(1), constants, 1);
            one_on_R_valid = true;
            if (isempty(out.xe)||out.ye(1)>0.01)
                one_on_R_valid = false;
            end
            counter = 1;
            crossed = false;
            alpha = sqrt(sigma)/(10*R);
            while one_on_R_valid
                % start by hunting downwards, until we cross
                if ~crossed
                    if counter>1
                        [~,out] = get_shape(psi_0_dot_iter_1_R(counter), constants, 1);
                    end
                    if isempty(out.xe)||out.ye(1)>0.01
                        % we've crossed, bracket
                        crossed = true;
                        psi_0_dot_lower_1_R = psi_0_dot_iter_1_R(counter);
                        psi_0_dot_upper_1_R = psi_0_dot_iter_1_R(counter-1);
                    else
                        alpha = alpha*2;
                        psi_0_dot_iter_1_R(counter+1) = psi_0_dot_iter_1_R(counter) - alpha; 
                    end
                else
                    %if we have crossed, now bisect
                    psi_0_dot_iter_1_R(counter) = (psi_0_dot_upper_1_R+psi_0_dot_lower_1_R)/2;
                    [~,out] = get_shape(psi_0_dot_iter_1_R(counter), constants, 1);
                    if isempty(out.xe)||out.ye(1)>0.01
                        % we're below the limit, current psi_0_dot is new upper
                        psi_0_dot_lower_1_R = psi_0_dot_iter_1_R(counter);
                    else
                        % we're above the limit, current psi_0_dot is new lower
                        psi_0_dot_upper_1_R = psi_0_dot_iter_1_R(counter);
                    end
                    alpha = abs(psi_0_dot_upper_1_R-psi_0_dot_lower_1_R);
                end
        
                s_end_iter_1_R(counter) = out.x(end);
                physical_solution_1_R(counter) = ~(isempty(out.xe)||out.ye(1)>0.01);
        
                if debug_on
                % for debugging
%                 fprintf('s iter = %0.2g\n',s_end_iter_1_R(counter))
%                 fprintf('psi 0 iter = %0.2g\n',psi_0_dot_iter_1_R(counter))
%                 fprintf('alpha = %0.2g\n',alpha)
                text(psi_0_dot_iter_1_R(counter), s_end_iter_1_R(counter), sprintf('%0.2g', counter))
                ylim([0,max(s_end_iter)])
                end
        
                if alpha<1e-8
                    break
                else
                    counter = counter + 1;
                end
            end
        end
    end
    
    
    % the final result is the largest s which is physical, unless we know that
    % we can't find a physical solution
    if greater_than_d_solution
        if consider_up_solution_greater_than_d
            psi_0_dot_end_up = outputs_up(1);
        end
        if consider_down_solution_greater_than_d
            psi_0_dot_end_down = outputs_down(1);
        end
    elseif can_find_solution
        [~, s_index] = max(s_end_iter(physical_solution));
        psi_0_dot_physical = psi_0_dot_iter(physical_solution);
        psi_0_dot_end = psi_0_dot_physical(s_index);
        
        if consider_up_solution
            [~, s_index] = max(s_end_iter(physical_solution_up));
            psi_0_dot_physical = psi_0_dot_iter_up(physical_solution_up);
            psi_0_dot_end_up = psi_0_dot_physical(s_index);
        end

        if one_on_R_valid
            [~, s_index] = max(s_end_iter_1_R(physical_solution_1_R));
            psi_0_dot_physical = psi_0_dot_iter_1_R(physical_solution_1_R);
            psi_0_dot_end_1_R = psi_0_dot_physical(s_index);
        end
    else
        psi_0_dot_end = psi_0_dot_iter(1);
    end
    
    % psi_0_dot_end = -1.02;
    if greater_than_d_solution
        % use condition 3 so that we end at r_max = d/2
        if consider_up_solution_greater_than_d && consider_down_solution_greater_than_d
            [~,out_down] = get_shape(psi_0_dot_end_down, constants, 3);
            [~,out_up] = get_shape(psi_0_dot_end_up, constants, 3);
            E_down = out_down.y(9,end);
            E_up = out_up.y(9,end);
            if E_down<E_up
                out = out_down;
            else
                out = out_up;
            end
        elseif consider_up_solution_greater_than_d
            [~,out] = get_shape(psi_0_dot_end_up, constants, 3);
        else
            [~,out] = get_shape(psi_0_dot_end_down, constants, 3);
        end
    elseif consider_up_solution
        [~,out_down] = get_shape(psi_0_dot_end, constants, 1);
        [~,out_up] = get_shape(psi_0_dot_end_up, constants, 1);
        E_down = out_down.y(9,end);
        E_up = out_up.y(9,end);
        E_1_R = abs(E_up*100);
        if one_on_R_valid
            [~,out_1_R] = get_shape(psi_0_dot_end_1_R, constants, 1);
            E_1_R = out_1_R.y(9,end);
        end
        [~, index] = min([E_down,E_up,E_1_R]);
    
        if index==1
            out = out_down;
        elseif index==2
            out = out_up;
        else
            out = out_1_R;
        end
    else
        [~,out] = get_shape(psi_0_dot_end, constants, 1);
    end
    % isempty(out.xe)
    
%     se = out.xe;
    % s = out.x;
    r = out.y(2,:);
    h = out.y(3,:);

%     E_down
%     E_up
    
    % solution = deval(out, linspace(0,out.x(end), 1000));
    % r = solution(2,:);
    % h = solution(3,:)-solution(3,end);
    
    if debug_on
        figure('Position',[400,100,800,600]);
        hold on
        axis equal
        xlabel('$r$')
        ylabel('$h$')
        plot(r, h-h(1), 'r-', 'displayname', 'free surface');
        t = linspace(-pi/2,pi/2,1000);
        x = cos(t)*R;
        % y = sin(t)*R+(R*cos(phi)+h(1));
        y = sin(t)*R+R*cos(phi);
        plot(x,y,'r:', 'displayname', 'microbead')
    end

end

%%
% 
function derivs = hamilton(s, y, const)
    psi = y(1);
    r = y(2);
    h = y(3);
    p_psi = y(4);
    p_r = y(5);
    p_h = y(6);
    E_bend = y(7);
    A = y(8);
    E_free = y(9);
    delA = y(10);

    lambda = const(1);

    f(1) = p_psi/(2*r)-sin(psi)/r;
    f(2) = cos(psi);
    f(3) = sin(psi);
    f(4) = cos(psi)*(p_psi/r-p_h)...
        +sin(psi)*p_r;
    f(5) = p_psi/r*(p_psi/(4*r)-sin(psi)/r)...
        +2/lambda^2;
    f(6) = 0;
    f(7) = p_psi^2/(4*r);
    f(8) = r;
    f(9) = f(7) + r*2/lambda^2*(1-cos(psi));
    f(10) = r.*(1-cos(psi));

    derivs = f';
end

% function derivs = hamilton(s, y, const)
%     psi = y(1);
%     r = y(2);
%     h = y(3);
%     p_psi = y(4);
%     p_r = y(5);
%     p_h = y(6);
%     E_bend = y(7);
%     A = y(8);
%     E_free = y(9);
% 
%     lambda = const(1);
% 
%     f(1) = p_psi/(2*r)-sin(psi)/r;
%     f(2) = cos(psi);
%     f(3) = sin(psi);
%     f(4) = cos(psi)*(p_psi/r-p_h)...
%         +sin(psi)*(2/lambda^2*r+p_r);
%     f(5) = p_psi/r*(p_psi/(4*r)-sin(psi)/r)...
%         +2/lambda^2*(1-cos(psi));
%     f(6) = 0;
%     f(7) = p_psi^2/(4*r);
%     f(8) = r;
%     f(9) = f(7) + r*2/lambda^2*(1-cos(psi));
% 
%     derivs = f';
% end

function [value,isterminal,direction] = myEventFcn(s,y,const)
    psi = y(1);
    r = y(2);
    h = y(3);
    p_psi = y(4);
    p_r = y(5);
    p_h = y(6);
    E_bend = y(7);
    A = y(8);
    E_free = y(9);
    delA = y(10);

    lambda = const(1);
    condition = const(2);
    d = const(3);

    if condition==1
        value = [psi, psi-pi];
        isterminal = [1,1];
        direction = [0,0];
    elseif condition==2
        value = [psi, psi-pi];
        isterminal = [1,1];
        direction = [0,0];
    elseif condition==3
        value = [psi, r-d/2];
        isterminal = [1,1];
        direction = [0,0];
    end

%     value = [y(1)];
%     isterminal = [1];
%     direction = [0];
end

function [diffs,out] = get_shape(psi_dot_init, constants, condition)

    lambda = constants(1);
    R = constants(2);
    phi = constants(3);
    d = constants(4);

    const(1) = lambda;
    const(2) = condition;
    const(3) = d;
    
    r_phi = R*sin(phi);
    psi_init = phi;
    p_psi_init = 2*r_phi*(psi_dot_init+sin(psi_init)/r_phi);
    p_r_init = -(p_psi_init^2/(4*r_phi)-p_psi_init*sin(psi_init)/r_phi...
        -2*r_phi/lambda^2)/cos(psi_init);
%     p_r_init = sin(phi)/(R*cos(phi))*(1+2*R^2/lambda^2-R^2*psi_dot_init^2);

    init_vals = [psi_init, r_phi, R*sin(phi), p_psi_init, p_r_init, 0,0,0,0,0];
    event_func = @(s,y) myEventFcn(s,y,const);
%     options = odeset('Events', @myEventFcn);
    options = odeset('Events', event_func, 'RelTol', 1e-7, 'AbsTol',1e-9);
    out = ode45(@(s, y) hamilton(s, y,const),[0,2*d],init_vals,...
        options);
    
    if condition==1
        psi_end = out.y(1,end);
        diffs = psi_end^2;
    elseif condition==2
        r_end = out.y(2,end);
        diffs = (r_end-d/2)^2;
    elseif condition==3
        r_end = out.y(2,end);
        diffs = r_end-d/2;
    elseif condition==4   
        psi_end = out.y(1,end);
        r_end = out.y(2,end);
        diffs = (r_end-d/2)^2+50*psi_end^2;
    end

%     if psi_end>pi/2
%         diffs = diffs*1000;
%     end
end