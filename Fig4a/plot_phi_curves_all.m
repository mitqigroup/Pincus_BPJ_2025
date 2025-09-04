function fig_handle = plot_phi_curves_all(data_structure, slice, phi_vals, varargin)

% parse inputs and set defaults
args = varargin;
nargs = numel(args);
k = 1;
fig_input = false;
plot_minima = false;
plot_local_minima = false;
use_newcolours = false;
remove_unphysical = false;
use_high_tension_approx = false;
use_aBaA_approx = false;
plot_relative = false;
ExtraCurveIndex = 0;
analytical_curves = 0;
plot_values = 1:6;
xaxis_name = '$\phi$';
yaxis_name = '$\Delta E$';
anno_string = [];
scaling = 1;
while k<=nargs
    % if we get a figure handle, use it
    if strcmpi(args{k},'FigureHandle')||strcmpi(args{k},'Figure_Handle')
        fig_handle_in = args{k+1};
        fig_input = true;
        k = k+1;
    elseif strcmpi(args{k},'xlabel')
        xaxis_name = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'ylabel')
        yaxis_name = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'AnalysisData')||strcmpi(args{k},'Analysis_Data')
        analysis_output = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'PlotMInima')||strcmpi(args{k},'Plot_minima')
        plot_minima = args{k+1};
        % we also want to make sure that we have analysis_output!
        if ~exist("analysis_output", 'var')
            error("to plot minima, we need the analysis output, make sure it's input first!")
        end
        k = k+1;
    elseif strcmpi(args{k},'PlotLocalMInima')||strcmpi(args{k},'Plot_local_minima')
        plot_local_minima = args{k+1};
        % we also want to make sure that we have analysis_output!
        if ~exist("analysis_output", 'var')
            error("to plot local minima, we need the analysis output, make sure it's input first!")
        end
        k = k+1;
    elseif strcmpi(args{k},'ColourOrder')||strcmpi(args{k},'ColorOrder')
        newcolors = args{k+1};
        use_newcolours = true;
        k = k+1;
    elseif strcmpi(args{k},'AnnotationString')
        anno_string = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'PlotRelative')
        plot_relative = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'RemoveUnphysical')
        remove_unphysical = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'ExtraCurveIndex')
        ExtraCurveIndex = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'AnalyticalCurves')
        analytical_curves = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'PlotCurves')
        plot_values = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'HighTensionApprox')||strcmpi(args{k},'HighTensionApproximation')
        use_high_tension_approx = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'aBaAApprox')||strcmpi(args{k},'aBaA_Approx')
        use_aBaA_approx = args{k+1};
        k = k+1;
    else
        error('intput not recognised!');
    end
    k = k+1;
end

ds = data_structure;

if fig_input
    fig_handle = figure(fig_handle_in,'Position',[400,100,800,600]);
else
    fig_handle = figure('Position',[400,100,800,600]);
    annotation('textbox', 'String', anno_string)
end
hold on
% xlim([0,180])
xlabel(xaxis_name)
ylabel(yaxis_name)
% lines = ["-", ":", ":", "--", ":", "--"];
colours = ["#016D24","#FE6100","#DC267F","#648FFF","#1D0C6B","#27E0D8","#FFB000"];
if use_newcolours
    colours = newcolors;
end
% jj = 16;
counter = 0;
jj = slice;
counter = counter+1;
epsilon = ds.parameter_set(jj,1);
n0 = ds.parameter_set(jj,2);
d = ds.parameter_set(jj,3);
R = ds.parameter_set(jj,4);
kD = ds.parameter_set(jj,5);
kappa = ds.parameter_set(jj,6);
alpha_i = ds.parameter_set(jj,7);
sigma = pi*R^2/d^2;
zeta = -epsilon*n0/kD;
E_unwrapped = [kD/2*alpha_i^2*d^2/(1+alpha_i),0,0,kD/2*alpha_i^2*d^2/(1+alpha_i),0,0,nan];

legend_names = {'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
            '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
            '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'};

E_total = ds.E_all_nonlinear(1,:,jj);
E_adhesion = ds.E_all_nonlinear(2,:,jj);
E_stretch_A = ds.E_all_nonlinear(3,:,jj);
E_stretch_B = ds.E_all_nonlinear(4,:,jj);
E_bend_A = ds.E_all_nonlinear(5,:,jj);
E_bend_B = ds.E_all_nonlinear(6,:,jj);

z_vals = (1-cos(phi_vals));
for ii=plot_values

    if ii==1
        primary = plot(rad2deg(phi_vals), ds.E_all_nonlinear(ii,:,jj)-E_unwrapped(ii),...
            '-','displayname', legend_names{ii}, 'Color', colours(ii));
        primary.LineWidth = 3;
    else
        p1 = plot(rad2deg(phi_vals), ds.E_all_nonlinear(ii,:,jj)-E_unwrapped(ii),...
            '-','displayname', legend_names{ii}, 'Color', colours(ii));
    end
end

switch ExtraCurveIndex
    case 1
        p2 = plot(rad2deg(phi_vals), ds.E_all_nonlinear(2,:,jj)+ ds.E_all_nonlinear(5,:,jj),...
            '-','displayname', '$E_\mathrm{adhesion}+E_\mathrm{bend,A}$', 'Color', colours(7));
    case 2
        p2 = plot(rad2deg(phi_vals), ds.E_all_nonlinear(2,:,jj)+ ds.E_all_nonlinear(4,:,jj),...
            '-','displayname', '$E_\mathrm{adhesion}+E_\mathrm{stretch,B}$', 'Color', colours(7));
    case 3
        p2 = plot(rad2deg(phi_vals),...
            ds.E_all_nonlinear(2,:,jj)+ ds.E_all_nonlinear(4,:,jj)...
            +ds.E_all_nonlinear(5,:,jj)-E_unwrapped(1),...
            '-','displayname', '$E_\mathrm{adhesion}+E_\mathrm{stretch,B}+E_\mathrm{bend,A}$',...
            'Color', colours(7));
    case 4
        xe = phi_vals;
        ye = 2*pi*R^2*epsilon*n0*(1-cos(xe'));
        pextra1 = plot(rad2deg(xe), ye, '-', 'MarkerSize',6,...
            'DisplayName','C1 $\equiv$ $2 \pi R^2 \epsilon n_0 (1-\cos \phi)$',...
            'Color', colours(2));
        ye = ye + 4*pi*kappa*(1-cos(xe'));
        pextra2 = plot(rad2deg(xe), ye, '-', 'MarkerSize',6,...
            'DisplayName','C2 $\equiv$ C1 $+ 4 \pi \kappa (1-\cos \phi)$',...
            'Color', colours(4));
        ab = (1+sigma*(2*(1-cos(xe'))-sin(xe').^2))*(1+alpha_i)-1;
        ye = ye + kD/2*(d^2-pi*R^2*sin(xe').^2).*ab.^2./(1+ab)-kD/2*d^2*alpha_i^2/(1+alpha_i);
        pextra3 = plot(rad2deg(xe), ye, '-', 'MarkerSize',6,...
            'DisplayName',"C3 $\equiv$ C2 $+$ annotation",...
            'Color', colours(5));
        annotation('textbox', 'String', ["$\frac{k_D}{2} \left(d^2 -\pi R^2 \sin^2 \phi \right) \frac{\alpha_B'^2}{1+\alpha_B'^2}$",...
            "$\alpha_B' = \left(1+ \sigma \left[ 2 \left(1-\cos \phi \right) - \sin^2 \phi \right] \right) (1+\alpha_i)-1$"])
        delA = ds.S_B_vals_nonlinear(:,jj)-(d^2-pi*R^2*sin(xe').^2);
%         ye = ye + ds.Sigma_vals_nonlinear(:,jj).*delA;
        ye = ye + kD*ab.*delA;
        pextra3 = plot(rad2deg(xe), ye, '-', 'MarkerSize',6,...
            'DisplayName',"C4 $\equiv$ C3 +$k_D \alpha_B' \Delta A$ ",...
            'Color', colours(7));
        ye = ye + ds.E_all_nonlinear(6,:,jj)';
        pextra3 = plot(rad2deg(xe), ye, '-', 'MarkerSize',6,...
            'DisplayName',"C5 $\equiv$ C4 +$E_\mathrm{bend,B}$ ",...
            'Color', colours(6));
    case 5
        p2 = plot(rad2deg(phi_vals),...
            E_stretch_A+E_stretch_B+E_bend_B,...
            '-','displayname', '$E_\mathrm{stretch,A}+E_\mathrm{stretch,B}+E_\mathrm{bend,B}$',...
            'Color', colours(7));
    case 6
        p2 = plot(rad2deg(phi_vals),...
            E_adhesion+E_bend_A,...
            '-','displayname', '$E_\mathrm{adhesion}+E_\mathrm{bend,A}$',...
            'Color', colours(7));
    case 7
        p2 = plot(rad2deg(phi_vals),...
            E_stretch_B+E_bend_B,...
            '-','displayname', '$E_\mathrm{stretch,B}+E_\mathrm{bend,B}$',...
            'Color', colours(7));

end

range_plot = 1:3:length(phi_vals); 
x = phi_vals(range_plot);
if any(analytical_curves==1)
    p3 = plot(rad2deg(x), 2*pi*R^2*epsilon*n0*(1-cos(x')), 'k^', 'MarkerSize',6,...
        'DisplayName','$2 \pi R^2 \epsilon n_0 (1-\cos \phi)$');
end
if any(analytical_curves==2)
    p4 = plot(rad2deg(x), 4*pi*kappa*(1-cos(x')), 'ko', 'MarkerSize',6,...
        'DisplayName','$4 \pi \kappa (1-\cos \phi)$');
end
if any(analytical_curves==3)
    ab = (1+sigma*(2*(1-cos(x))-sin(x).^2))*(1+alpha_i)-1;
    ytest = kD/2*(d^2-pi*R^2*sin(x).^2).*ab.^2./(1+ab)-kD/2*d^2*alpha_i^2/(1+alpha_i);
    p5 = plot(rad2deg(x), ytest, 'ks', 'MarkerSize',6,...
        'DisplayName','See above');
    annotation('textbox', 'String', ["$\frac{k_D}{2} \left(d^2 -\pi R^2 \sin^2 \phi \right) \frac{\alpha_B'^2}{1+\alpha_B'^2}$",...
        "$\alpha_B' = \left(1+ \sigma \left[ 2 \left(1-\cos \phi \right) - \sin^2 \phi \right] \right) (1+\alpha_i)-1$"])
end
if any(analytical_curves==4)
    % get delta A by subtracting away spherical cap area
    ab = (1+sigma*(2*(1-cos(x))-sin(x).^2))*(1+alpha_i)-1;
    delA = ds.S_B_vals_nonlinear(range_plot,jj)'-(d^2-pi*R^2*sin(x).^2);
    ytest = kD*ab.*delA;
    p5 = plot(rad2deg(x), ytest, 'k>', 'MarkerSize',6,...
        'DisplayName',"$k_D \alpha_B' \Delta A$");
end
if any(analytical_curves==5)
    range_plot = 1:length(phi_vals); 
    x = phi_vals(range_plot);
    ab = (1+sigma*(2*(1-cos(x))-sin(x).^2))*(1+alpha_i)-1;
    ytest1 = kD/2*(d^2-pi*R^2*sin(x).^2).*ab.^2./(1+ab)-kD/2*d^2*alpha_i^2/(1+alpha_i);
    p5 = plot(rad2deg(x), ytest1, '--', 'MarkerSize',6,...
        'DisplayName',"$k_D \alpha_B' (S_\mathrm{B}-\Delta S_\mathrm{B})$", 'Color',colours(4));
    annotation('textbox', 'String',...
        "$\alpha_B' = \left(1+ \sigma \left[ 2 \left(1-\cos \phi \right) - \sin^2 \phi \right] \right) (1+\alpha_i)-1$")
    % get delta A by subtracting away spherical cap area
    ab = (1+sigma*(2*(1-cos(x))-sin(x).^2))*(1+alpha_i)-1;
    delA = ds.S_B_vals_nonlinear(range_plot,jj)'-(d^2-pi*R^2*sin(x).^2);
    ytest2 = kD*ab.*delA;
    p5 = plot(rad2deg(x), ytest2, ':', 'MarkerSize',6,...
        'DisplayName',"$k_D \alpha_B' \Delta S_\mathrm{B}$", 'Color',colours(4));
    % finally plot all others
    p2 = plot(rad2deg(phi_vals), ds.E_all_nonlinear(2,:,jj) + ...
        ds.E_all_nonlinear(5,:,jj) + ds.E_all_nonlinear(3,:,jj)+...
        ds.E_all_nonlinear(4,:,jj)-ytest1-ytest2-E_unwrapped(1),...
        '-','displayname', 'All other $E$', 'Color', colours(7));
end

if use_high_tension_approx
    constants = [epsilon, n0, d, R, kD, kappa, alpha_i];
    [delE, ~] = high_tension_approx(phi_vals, constants);

    plot(rad2deg(phi_vals), delE, ':', 'Color', [0.5,0.5,0.5],...
        'DisplayName','High Tension Approximation')
end

if use_aBaA_approx
    constants = [epsilon, n0, d, R, kD, kappa, alpha_i];
    [delE, ~] = aBaA_approx(phi_vals, constants);

    plot(rad2deg(phi_vals), delE, '--', 'Color', [0.5,0.5,0.5],...
        'DisplayName','$\alpha_A = \alpha_B$, no free membrane')
end

if plot_relative
    lines = fig_handle.Children(1).Children;
    for ii=1:length(lines)
        y = lines(ii).YData;
        ymax = max(abs(y));
        y = y/ymax;
        fig_handle.Children(1).Children(ii).YData = y;
        if ii==6
            scaling = ymax;
        end
    end
    ylabel('$\Delta E/\max(|\Delta E|)$')
end

if plot_local_minima
%     scatter(rad2deg(analysis_output.phi_at_min), ...
%         (analysis_output.E_all_min(1,:)-E_unwrapped(:,1)')./scaling,50,...
%         'p','MarkerEdgeColor',p1.Color, 'MarkerFaceColor',p1.Color,...
%         'HandleVisibility','off')
    scatter(rad2deg(analysis_output.phi_loc_min), ...
        (analysis_output.E_loc_min(1,:)-E_unwrapped(:,1)')./scaling,50,...
        'o','MarkerEdgeColor',primary.Color, 'MarkerFaceColor',primary.Color,...
        'HandleVisibility','off')
    % also put dots at the endpoints if they're minima
    counter = 0;
    for jj = slice
        counter = counter+1;
        if analysis_output.zero_phi_is_min(counter)
            scatter(rad2deg(min(phi_vals)), ...
                (analysis_output.E_low_phi(counter)-E_unwrapped(counter,1)')./scaling(counter),50,...
                'o','MarkerEdgeColor',primary.Color, 'MarkerFaceColor',primary.Color,...
                'HandleVisibility','off')
        end
        if analysis_output.pi_phi_is_min(counter)
            scatter(rad2deg(max(phi_vals)), ...
                (analysis_output.E_high_phi(counter)-E_unwrapped(counter,1)')./scaling(counter),50,...
                'o','MarkerEdgeColor',primary.Color, 'MarkerFaceColor',primary.Color,...
                'HandleVisibility','off')
        end
    end
end

if plot_minima
    scatter(rad2deg(analysis_output.phi_at_min(counter)), ...
        (analysis_output.E_all_min(1,counter)-E_unwrapped(1))./scaling,300,...
        'p','MarkerEdgeColor',primary.Color, 'MarkerFaceColor',primary.Color,...
        'HandleVisibility','off')
end

xlim([0,180])
legend;

end