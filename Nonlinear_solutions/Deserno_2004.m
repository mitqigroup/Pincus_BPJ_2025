clear variables
% warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
warning('on', 'MATLAB:ode45:IntegrationTolNotMet');
R = 1;
% z = eps(1);
% sigma = 0.7004;
% sigma = 0.399999975000000;
% sigma = 5;
% w = 4+sigma*2+0.05;
% w = 1.2125*4;
alpha_i = 1e-05;
kD = 0.3;
% epsilon = 0.00052*kD;
epsilon = 8.5441e-05;
kappa = 1e-19*1e12;     % picoJ
w = epsilon*2*R^2/kappa;
% w = 5;
% w = 1e-4*2*R^2/kappa;
% sigma = 0.00050782*R^2/kappa;
% sigma = 0.0005*R^2/kappa;
% sigma = kD*alpha_i*R^2/kappa;
% sigma = 5.340218205595971e-04*R^2/kappa;
% sigma = 12.221964;
sigma = 1000;
% sigma_vals = [5.9979    6.0306    6.0634    6.0947    6.1252,...
%     6.1548    6.1837    6.2116    6.2389    6.2653    6.2909,...
%     6.3157    6.3398    6.3627    6.3861    6.4123    6.4123,...
%     6.4112    6.4111    6.4106    6.4102    6.4094    6.4085,...
%     6.4075    6.4063    6.4054    6.4043    6.4030,...
%     6.4017    6.4005];
% sigma = 5000;
% lambda = R/sqrt(sigma);

debug_on = false;
% debug_on = true;

% phi_vals = flip(deg2rad(linspace(130, 179,10)));
% phi_vals = deg2rad(linspace(140,160,30));
% phi_vals = deg2rad(linspace(148.5, 152,10));
% phi_vals = deg2rad(linspace(148.6, 148.8,10));
phi_vals = deg2rad(linspace(0.01, 179.9,100));
% phi_vals = 2.814430228410981;
% phi_vals = deg2rad(149);
% phi_vals = deg2rad(30);
% phi_vals = deg2rad([5,50,80])
z_vals = 1-cos(phi_vals);
% z_vals = linspace(1.84,1.87,30);
% z_vals = linspace(0+eps(1),1.9999, 50);
% z_vals = [1.949585884727274,1.959049245604295,1.967548105003720]


rend = zeros(size(z_vals));
rphi = zeros(size(z_vals));
A = zeros(size(z_vals));
E_free =  zeros(size(z_vals));
E_bend =  zeros(size(z_vals));

tic
for jj = 1:length(z_vals)
z = z_vals(jj)
phi = acos(1-z);
% sigma = sigma_vals(jj);
lambda = R/sqrt(sigma);
constants = [lambda, R, phi];

if debug_on
    % psi_0 = linspace(-sqrt(sigma)/R-1,1/R,200);
    psi_0 = linspace(-2*sqrt(sigma)/R-2/R,2/R,200);
    % psi_0 = linspace(-1.16,-1.13,500);
    % psi_0 = linspace(-1e6,1,500);
    s_end = zeros(size(psi_0));
    correct = zeros(size(psi_0));
    correct2 = zeros(size(psi_0));
    for ii = 1:length(psi_0)
        [~,out] = get_shape(psi_0(ii), constants, 1);
        s_end(ii) = out.x(end);
        correct(ii) = isempty(out.xe)||out.ye(1)>0.01;
        correct2(ii) = isempty(out.xe);
    end

    figure('Position',[400,100,800,600]);
    hold on
    plot(psi_0, s_end, '-')
    xlabel('$\dot{\psi_0}$')
    ylabel('$s_\mathrm{max}$')
    yyaxis right
    ylim([-0.05,1.05])
    plot(psi_0, correct)
    ylabel('$\phi =\pi$ True/False')
    annotation("textbox","String",{sprintf('$\\phi = %0.4g^{\\circ}$', rad2deg(phi)),...
        sprintf('$R = %0.2g$', R), sprintf('$\\bar{\\Sigma} = %0.2g$', sigma)});
    % xlim([-20*sqrt(sigma)/R-1,1/R])
    yyaxis left
end

psi_0_upper = 1/R;
psi_0_lower = -4*sqrt(sigma)/R-2/R;

N = 1000;
psi_0_iter = zeros(1,N);
s_end_iter = zeros(1,N);
physical_solution = false(1,N);
can_find_solution = true;

psi_0_iter_up = zeros(1,N);
s_end_iter_up = zeros(1,N);
physical_solution_up = false(1,N);

consider_up_solution = false;
one_on_R_valid = false;

if z<1
    psi_0_iter(1) = psi_0_lower;
    alpha = sqrt(sigma)/(10*R);
    counter = 1;
    crossed = false;
    while true
        % start by hunting upwards, until we cross
        if ~crossed
            [~,out] = get_shape(psi_0_iter(counter), constants, 1);
            if isempty(out.xe)||out.ye(1)>0.01
                % we've crossed, bracket
                crossed = true;
                psi_0_upper = psi_0_iter(counter);
                psi_0_lower = psi_0_iter(counter-1);
            else
                alpha = alpha*2;
                psi_0_iter(counter+1) = psi_0_iter(counter) + alpha; 
                if psi_0_iter(counter+1)>1/R
                    % we've hit the unphysical limit, bracket
                    crossed = true;
                    psi_0_upper = 1/R;
                    psi_0_lower = psi_0_iter(counter);
                end
            end
        else
            %if we have crossed, now bisect
            psi_0_iter(counter) = (psi_0_upper+psi_0_lower)/2;
            [~,out] = get_shape(psi_0_iter(counter), constants, 1);
            if isempty(out.xe)||out.ye(1)>0.01
                % we're above the limit, current psi_0 is new upper
                psi_0_upper = psi_0_iter(counter);
            else
                % we're below the limit, current psi_0 is new lower
                psi_0_lower = psi_0_iter(counter);
            end
            alpha = abs(psi_0_upper-psi_0_lower);
        end

        s_end_iter(counter) = out.x(end);
        physical_solution(counter) = ~(isempty(out.xe)||out.ye(1)>0.01);

        if debug_on
            % for debugging
            fprintf('s iter = %0.2g\n',s_end_iter(counter))
            fprintf('psi 0 iter = %0.2g\n',psi_0_iter(counter))
            fprintf('alpha = %0.2g\n',alpha)
            text(psi_0_iter(counter), s_end_iter(counter), sprintf('%0.2g', counter))
            ylim([0,max(s_end_iter)])
        end

        if alpha<1e-8
            break
        else
            counter = counter + 1;
        end
    end
else
%     psi_0_iter(1) = 1/(R);
%     psi_0_iter(1) = psi_0(find(correct==0,1,'last'));
    psi_0_iter(1) = -3/(10*R);
    % check that this first point is valid
    [~,out] = get_shape(psi_0_iter(1), constants, 1);
    times_searched = 0;
    counter = 1;
    crossed = false;
    alpha = sqrt(sigma)/(10*R);
    while (isempty(out.xe)||out.ye(1)>0.01)
        % if not valid, start at 1/R and search down
        psi_0_iter(1) = 1/R - alpha*times_searched/2;
        [~,out] = get_shape(psi_0_iter(1), constants, 1);
        times_searched = times_searched + 1;
        if psi_0_iter(1)<(-4*sqrt(sigma)/R-2/R)
            % at this point, just give up. I've got no idea what to do
            % here, maybe we can pick a close point in the phi-space? I
            % think for now just calculate at psi_0 = 1/R
            can_find_solution = false;
            psi_0_iter(1) = 1/R;
            [~,out] = get_shape(psi_0_iter(1), constants, 1);
            s_end_iter(counter) = out.x(end);
            physical_solution(counter) = ~(isempty(out.xe)||out.ye(1)>0.01);
            break
        end
    end
    while can_find_solution
        % start by hunting downwards, until we cross
        if ~crossed
            if counter>1
                [~,out] = get_shape(psi_0_iter(counter), constants, 1);
            end
            if isempty(out.xe)||out.ye(1)>0.01
                % we've crossed, bracket
                crossed = true;
                psi_0_lower = psi_0_iter(counter);
                psi_0_upper = psi_0_iter(counter-1);
            else
                alpha = alpha*2;
                psi_0_iter(counter+1) = psi_0_iter(counter) - alpha; 
            end
        else
            %if we have crossed, now bisect
            psi_0_iter(counter) = (psi_0_upper+psi_0_lower)/2;
            [~,out] = get_shape(psi_0_iter(counter), constants, 1);
            if isempty(out.xe)||out.ye(1)>0.01
                % we're below the limit, current psi_0 is new upper
                psi_0_lower = psi_0_iter(counter);
            else
                % we're above the limit, current psi_0 is new lower
                psi_0_upper = psi_0_iter(counter);
            end
            alpha = abs(psi_0_upper-psi_0_lower);
        end

        s_end_iter(counter) = out.x(end);
        physical_solution(counter) = ~(isempty(out.xe)||out.ye(1)>0.01);

        if debug_on
        % for debugging
        fprintf('s iter = %0.2g\n',s_end_iter(counter))
        fprintf('psi 0 iter = %0.2g\n',psi_0_iter(counter))
        fprintf('alpha = %0.2g\n',alpha)
        text(psi_0_iter(counter), s_end_iter(counter), sprintf('%0.2g', counter))
        ylim([0,max(s_end_iter)])
        end

        if alpha<1e-8
            break
        else
            counter = counter + 1;
        end
    end
    % we also want to hunt up, just to be sure
    counter = 1;
    crossed = false;
    psi_0_iter_up(1) = psi_0_iter(1);
    alpha = sqrt(sigma)/(30*R);
    while can_find_solution
        % start by hunting upwards, until we cross or get past 1/R
        if ~crossed
            [~,out] = get_shape(psi_0_iter_up(counter), constants, 1);
            if psi_0_iter_up(counter)>1/R
                % check that the 1/R solution is feasible, if it is, we
                % should stop the search, if it isn't, we should bisect!
                [~,out_test] = get_shape(1/R, constants, 1);
                if isempty(out_test.xe)||out_test.ye(1)>0.01
                    % move to bisection
                    psi_0_upper = 1/R;
                    psi_0_lower = psi_0_iter_up(counter-1);
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
                psi_0_upper = psi_0_iter_up(counter);
                psi_0_lower = psi_0_iter_up(counter-1);
            else
%                 alpha = alpha*2;
                psi_0_iter_up(counter+1) = psi_0_iter_up(counter) + alpha; 
            end
        else
            %if we have crossed, now bisect
            psi_0_iter_up(counter) = (psi_0_upper+psi_0_lower)/2;
            [~,out] = get_shape(psi_0_iter_up(counter), constants, 1);
            if isempty(out.xe)||out.ye(1)>0.01||(out.ye(1)>0.01&&out.ye(2)>d/2)
                % we're below the limit, current psi_0 is new upper
                psi_0_upper = psi_0_iter_up(counter);
            else
                % we're above the limit, current psi_0 is new lower
                psi_0_lower = psi_0_iter_up(counter);
            end
            alpha = abs(psi_0_upper-psi_0_lower);
        end

        s_end_iter_up(counter) = out.x(end);
        physical_solution_up(counter) = ~(isempty(out.xe)||out.ye(1)>0.01);

        if debug_on
            % for debugging
            fprintf('s iter = %0.2g\n',s_end_iter_up(counter))
            fprintf('psi 0 iter = %0.2g\n',psi_0_iter_up(counter))
            fprintf('alpha = %0.2g\n',alpha)
            text(psi_0_iter_up(counter), s_end_iter_up(counter), sprintf('%0.2g', counter), 'color', 'red')
%             ylim([0,max(s_end_iter)])
        end

        if alpha<1e-8
            consider_up_solution = true;
            break
        else
            counter = counter + 1;
        end
    end
    % finally, we always want to check the 1/R point if we haven't already,
    % then hunt downwards
    if psi_0_iter(1) ~= 1/R
        psi_0_iter_1_R(1) = 1/R;
        % check that this first point is valid
        [~,out] = get_shape(psi_0_iter_1_R(1), constants, 1);
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
                    [~,out] = get_shape(psi_0_iter_1_R(counter), constants, 1);
                end
                if isempty(out.xe)||out.ye(1)>0.01
                    % we've crossed, bracket
                    crossed = true;
                    psi_0_lower_1_R = psi_0_iter_1_R(counter);
                    psi_0_upper_1_R = psi_0_iter_1_R(counter-1);
                else
                    alpha = alpha*2;
                    psi_0_iter_1_R(counter+1) = psi_0_iter_1_R(counter) - alpha; 
                end
            else
                %if we have crossed, now bisect
                psi_0_iter_1_R(counter) = (psi_0_upper_1_R+psi_0_lower_1_R)/2;
                [~,out] = get_shape(psi_0_iter_1_R(counter), constants, 1);
                if isempty(out.xe)||out.ye(1)>0.01
                    % we're below the limit, current psi_0 is new upper
                    psi_0_lower_1_R = psi_0_iter_1_R(counter);
                else
                    % we're above the limit, current psi_0 is new lower
                    psi_0_upper_1_R = psi_0_iter_1_R(counter);
                end
                alpha = abs(psi_0_upper_1_R-psi_0_lower_1_R);
            end
    
            s_end_iter_1_R(counter) = out.x(end);
            physical_solution_1_R(counter) = ~(isempty(out.xe)||out.ye(1)>0.01);
    
            if debug_on
            % for debugging
            fprintf('s iter = %0.2g\n',s_end_iter_1_R(counter))
            fprintf('psi 0 iter = %0.2g\n',psi_0_iter_1_R(counter))
            fprintf('alpha = %0.2g\n',alpha)
            text(psi_0_iter_1_R(counter), s_end_iter_1_R(counter), sprintf('%0.2g', counter))
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
if can_find_solution
    [~, s_index] = max(s_end_iter(physical_solution));
    psi_0_physical = psi_0_iter(physical_solution);
    psi_0_end = psi_0_physical(s_index);
    
    if consider_up_solution
        [~, s_index] = max(s_end_iter(physical_solution_up));
        psi_0_physical = psi_0_iter_up(physical_solution_up);
        psi_0_end_up = psi_0_physical(s_index);
    end

    if one_on_R_valid
        [~, s_index] = max(s_end_iter_1_R(physical_solution_1_R));
        psi_0_physical = psi_0_iter_1_R(physical_solution_1_R);
        psi_0_end_1_R = psi_0_physical(s_index);
    end
else
    psi_0_end = psi_0_iter(1);
end
    
% figure('Position',[400,100,800,600]);
% hold on
% % plot(psi_0_iter, s_end_iter(2:end), ':')
% scatter(psi_0_iter, s_end_iter, [], 1:length(psi_0_iter))
% textscatter(psi_0_iter, s_end_iter, compose('%d', 1:length(psi_0_iter)),...
%     "TextDensityPercentage",50)
% plot(psi_0, s_end, '-')
% xlabel('$\psi_0$')
% ylabel('$s_\mathrm{max}$')
% yyaxis right
% ylim([-0.05,1.05])
% plot(psi_0, correct)
% colorbar
% annotation("textbox","String",{sprintf('$\\phi = %0.3g^{\\circ}$', rad2deg(phi)),...
%     sprintf('$R = %0.2g$', R), sprintf('$\\bar{\\sigma} = %0.2g$', sigma)});
% xlim([-4*sqrt(sigma)/R-2,1/R])
% yyaxis left

% psi_0_end = -1.02;
if consider_up_solution
    [~,out_down] = get_shape(psi_0_end, constants, 1);
    [~,out_up] = get_shape(psi_0_end_up, constants, 1);
    E_down = out_down.ye(9);
    E_up = out_up.ye(9);
    E_1_R = E_up*100;
    if one_on_R_valid
        [~,out_1_R] = get_shape(psi_0_end_1_R, constants, 1);
        E_1_R = out_1_R.ye(9);
    end
    [val, index] = min([E_down,E_up,E_1_R]);

    if index==1
        out = out_down;
    elseif index==2
        out = out_up;
    else
        out = out_1_R;
    end
%     if E_down<E_up
%         out = out_down;
%     else
%         out = out_up;
%     end
else
    [~,out] = get_shape(psi_0_end, constants, 1);
    index = 4;
end
which_solution(jj) = index;
% isempty(out.xe)

se = out.xe;
% s = out.x;
r = out.y(2,:);
h = out.y(3,:);
rend(jj) = r(end);
rphi(jj) = r(1);
A(jj) = out.y(8,end);
E_free(jj) = out.y(9,end);
E_bend(jj) = out.y(7,end);

% solution = deval(out, linspace(0,out.x(end), 1000));
% r = solution(2,:);
% h = solution(3,:)-solution(3,end);

if debug_on
    figure('Position',[400,100,800,600]);
    hold on
    axis equal
    xlabel('$r$')
    ylabel('$h$')
    t = linspace(-pi/2,pi/2,1000);
    x = cos(t)*R;
    % y = sin(t)*R+(R*cos(phi)+h(1));
    y = sin(t)*R+R*cos(phi);
    plot(x,y,'r:', 'displayname', 'microbead')

    if consider_up_solution
        r_down = out_down.y(2,:);
        h_down = out_down.y(3,:);
        plot(r_down, h_down-h_down(1), 'r-', 'displayname',...
            sprintf('down, $E = %0.5g$', E_down));

        r_up = out_up.y(2,:);
        h_up = out_up.y(3,:);
        plot(r_up, h_up-h_up(1), 'b-', 'displayname',...
            sprintf('up, $E = %0.5g$', E_up));
        if one_on_R_valid
            r_1_R = out_1_R.y(2,:);
            h_1_R = out_1_R.y(3,:);
            plot(r_1_R, h_1_R-h_1_R(1), 'g-', 'displayname',...
            sprintf('1 on R, $E = %0.5g$', E_1_R));
        end
    else
        plot(r, h-h(1), 'r-', 'displayname', 'free surface');
    end
    legend
end
end
toc

%%
% w = 5;

% 
delA = 2*pi*A - pi*(rend.^2-rphi.^2);
delA/R^2;
rend/R;
% 
% test1 = delA*2/lambda^2;
% test2 = (E_free-E_bend)*2*pi;
% % test3 = trapz(s, r.*(1-cos(out.y(1,:))))*2*pi*2/lambda^2
% 
% figure();
% hold on
% plot(z_vals, test1)
% plot(z_vals, test2)
% % plot(z_vals, test3)
% 
% 
E = -(w-4)*z_vals+sigma*z_vals.^2+E_free;
% 
E_des(1,:) = E;
E_des(2,:) = -w*z_vals;
E_des(3,:) = sigma*z_vals.^2;
E_des(4,:) = (E_free-E_bend);
E_des(5,:) = 4*z_vals;
E_des(6,:) = E_bend;
E_des(7,:) = E_des(3,:)+E_des(4,:);
E_des(8,:) = E_des(4,:)+E_des(6,:);
E_des = E_des*pi*kappa;
% 
% k_vals = sqrt(z_vals.*(2-z_vals));
% E_free_linear = R/lambda*(k_vals.^3./(1-k_vals.^2)).*besselk(0,k_vals*R/lambda)/besselk(1,k_vals*R/lambda);
% E_linear = E_free_linear - w*z_vals + sigma*z_vals.^2 + 4*z_vals;
% 
% figure();
% hold on
% xlabel('$\phi$')
% ylabel('$\bar{E}$')
% % plot(z_vals, E)
% % plot(z_vals, E_free)
% % plot(rad2deg(acos(1-z_vals)), E*pi*kappa)
% plot(rad2deg(acos(1-z_vals)), E/max(abs(E)))
% % plot(rad2deg(acos(1-z_vals)), sum(E_des(2:end-1,:)))
% % plot(rad2deg(acos(1-z_vals)), E_linear*pi*kappa)
% 
% figure();
% plot(z_vals, E)
% plot(rad2deg(acos(1-z_vals)), E_des)
% legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
%     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
%     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$', '$E_\mathrm{bend}$'}, 'Box','off',...
%     'location', 'best')

E_bend_high_ten = pi*kappa*R*sin(phi_vals)/lambda.*(4-2*sqrt(2-2*cos(phi_vals)).*cot(phi_vals/2));
E_stretch_high_ten = 4*pi*kappa*R*sin(phi_vals).*(2/lambda*(1-sqrt((1+cos(phi_vals))/2)));

figure();
hold on
plot(rad2deg(acos(1-z_vals)), E_des([4,6,8],:))
plot(rad2deg(acos(1-z_vals)), E_bend_high_ten)
plot(rad2deg(acos(1-z_vals)), E_stretch_high_ten)
legend({'$E_\mathrm{stretch,B}$','$E_\mathrm{bend,B}$','$E_\mathrm{total,B}$',...
    '$E_\mathrm{bend,approx}$', '$E_\mathrm{total,approx}$'}, 'Box','off',...
    'location', 'best')

% 
% % figure();
% % plot(rad2deg(acos(1-z_vals)), delA)
% 
% % % get minimum
% % [min_E, E_index] = min(E_des(1,:));
% % E_des_min = E_des(:,E_index);
% % figure();
% % x = [1e-4,1e-3];
% % y = repmat(E_des_min,1,2);
% % plot(x,y);
% % legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
% %     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
% %     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
% %     'location', 'best')

%% test single point
% phi = deg2rad(60);
constants = [lambda, R, phi];
% [diffs,out] = get_shape(-1.1, constants, 1);
[diffs,out] = get_shape(psi_0_end, constants, 1);
% isempty(out.xe)

% s = out.x;
% r = out.y(2,:);
% h = out.y(3,:);

solution = deval(out, linspace(0,out.x(end), 1000));
r = solution(2,:);
h = solution(3,:)-solution(3,end);

% figure('Position',[400,100,800,600]);
% hold on
% axis equal
% xlabel('$r$')
% ylabel('$h$')
% plot(r, h-h(1), 'r-', 'displayname', 'free surface');
% t = linspace(-pi/2,pi/2,1000);
% x = cos(t)*R;
% % y = sin(t)*R+(R*cos(phi)+h(1));
% y = sin(t)*R+R*cos(phi);
% plot(x,y,'r:', 'displayname', 'microbead')

%%
[diffs,out] = get_shape(psi_0_end, constants, 1);
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

    value = [y(1), y(1)-pi];
    isterminal = [1,1];
    direction = [0,0];

%     value = [y(1)];
%     isterminal = [1];
%     direction = [0];
end

function [diffs,out] = get_shape(psi_dot_init, constants, condition)

    lambda = constants(1);
    R = constants(2);
    phi = constants(3);

    const(1) = lambda;
    
    r_phi = R*sin(phi);
    psi_init = phi;
    p_psi_init = 2*r_phi*(psi_dot_init+sin(psi_init)/r_phi);
    p_r_init = -(p_psi_init^2/(4*r_phi)-p_psi_init*sin(psi_init)/r_phi...
        -2*r_phi/lambda^2)/cos(psi_init);

    init_vals = [psi_init, r_phi, R*sin(phi), p_psi_init, p_r_init, 0,0,0,0];
    event_func = @(s,y) myEventFcn(s,y,const);
%     options = odeset('Events', @myEventFcn);
    options = odeset('Events', event_func, 'RelTol', 1e-7, 'AbsTol',1e-9);
    out = ode45(@(s, y) hamilton(s, y,const),[0,1000],init_vals,...
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