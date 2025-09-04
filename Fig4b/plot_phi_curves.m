function fig_handle = plot_phi_curves(data_structure, slice, phi_vals, varargin)

% parse inputs and set defaults, starting at 4
args = varargin;
nargs = numel(args);
k = 1;
fig_input = false;
plot_minima = false;
use_newcolours = false;
remove_unphysical = false;
plot_local_minima = false;
include_legend = true;
use_rel_E = true;
xaxis_name = '$\phi$';
yaxis_name = '$\Delta E$';
anno_string = [];
ThickerLines = [];
scaling = ones(size(slice));
disp_index = 0;
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
    elseif strcmpi(args{k},'Relative_Energy')||strcmpi(args{k},'RelativeEnergy')
        use_rel_E = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'DisplayNameIndex')
        disp_index = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'AnnotationString')
        anno_string = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'IncludeLegend')||strcmpi(args{k},'Include_Legend')
        include_legend = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'RemoveUnphysical')
        remove_unphysical = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'ThickerLines')
        ThickerLines = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'ScaleFactor')||strcmpi(args{k},'Scale_Factor')
        scaling = args{k+1};
        k = k+1;
    else
        error('intput not recognised!');
    end
    k = k+1;
end

ds = data_structure;

outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds,slice,phi_vals,...
    'RemoveUnphysical', true, 'cleanData', true);
% if remove_unphysical
%     % we will need this analysis later
%     outputs_analyse = analyse_phi_curve(ds.E_all_nonlinear,ds,slice,phi_vals,...
%         'RemoveUnphysical', true);
% end

if fig_input
    fig_handle = figure(fig_handle_in);
else
    fig_handle = figure('Position',[400,100,800,600]);
    annotation('textbox', 'String', anno_string)
end
hold on
% xlim([0,180])
xlabel(xaxis_name)
ylabel(yaxis_name)
% lines = ["-", ":", ":", "--", ":", "--"];
% colours = ['k', 'b', 'r', 'r', 'g', 'g'];
if use_newcolours
    colororder(newcolors)
end
% jj = 16;
counter = 0;
for jj = slice
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
    E_unwrapped(counter,:) = [kD/2*alpha_i^2*d^2/(1+alpha_i),0,0,kD/2*alpha_i^2*d^2/(1+alpha_i),0,0,nan];
    
    switch disp_index
        case 1
            display_name = sprintf('$%s, \\epsilon = %0.3g$', ...
                get_curve_name(outputs_analyse.curve_type(counter)), epsilon);
        case 3
            display_name = sprintf('$%s, d = %0.3g$', ...
                get_curve_name(outputs_analyse.curve_type(counter)), d);
        case 4
            display_name = sprintf('$%s, R = %0.3g$', ...
                get_curve_name(outputs_analyse.curve_type(counter)), R);
        case 5
            display_name = sprintf('$%s, k_D = %0.3g$', ...
                get_curve_name(outputs_analyse.curve_type(counter)), kD);
        case 6
            display_name = sprintf('$%s, \\kappa = %0.3g$', ...
                get_curve_name(outputs_analyse.curve_type(counter)), kappa);
        case 7
            display_name = sprintf('$%s, \\alpha_i = %0.3g$', ...
                get_curve_name(outputs_analyse.curve_type(counter)), alpha_i);
        case 8
            display_name = sprintf('$%s, \\sigma = %0.3g$', ...
                get_curve_name(outputs_analyse.curve_type(counter)), sigma);
        case 9
            display_name = sprintf('$%s, \\zeta = %0.3g$', ...
                get_curve_name(outputs_analyse.curve_type(counter)), zeta);
        case 10
            display_name = sprintf('$%s, \\bar{\\epsilon} = %0.3g$', ...
                get_curve_name(outputs_analyse.curve_type(counter)), 2*epsilon*n0*R^2/kappa);
        case 11
            display_name = sprintf('Type %i', outputs_analyse.curve_type(counter));
        case 12
            display_name = sprintf('$%s$', get_curve_name(outputs_analyse.curve_type(counter)));
        otherwise
            display_name = ' ';
    end
    
    ab = ds.alpha_B_vals_nonlinear(:,jj);
    E_n = ds.E_all_nonlinear(1,:,jj);
    if remove_unphysical
        unphysical_points = outputs_analyse.unphysical_points{counter};
        E_n(unphysical_points) = nan;
        % get these points and those on either side, and add them to
        % another cell array to plot at the end
    end
    %     if 1==1||1==4||1==6
    %         E_t = E_all_toroid(1,:,jj);
    %         E_n(E_n>E_t) = nan;
    %     end
    E_n = E_n - E_unwrapped(counter,1);
    if use_rel_E
        scaling(counter) = max(abs(E_n));
    end
    E_n = E_n/scaling(counter);
    limit_alpha_B = ab>=0.1;
    too_large_alpha_B = ab>0.04;
    E_ab = E_n;
    E_ab(~too_large_alpha_B) = nan;
    E_ab(limit_alpha_B) = nan;
    E_n(too_large_alpha_B) = nan;
    dont_plot_max(counter) = any(limit_alpha_B);
    if use_newcolours
        p1 = plot(rad2deg(phi_vals), E_n,...
            '-','displayname', display_name, 'Color',newcolors(counter));
        plot(rad2deg(phi_vals), E_ab,...
            '--','displayname', display_name, 'Color',newcolors(counter), 'HandleVisibility','off');
    else
        p1 = plot(rad2deg(phi_vals), E_n,...
            '-','displayname', display_name);
        plot(rad2deg(phi_vals), E_ab,...
            '--','displayname', display_name, 'HandleVisibility','off');
    end
    
    if any(ThickerLines==counter)
        p1.LineWidth = 3;
    end

    if remove_unphysical&&~isempty(unphysical_points)
        % get points near unphysical
        unphys_up = unphysical_points+1;
        unphys_up = unphys_up(unphys_up<=length(phi_vals));
        unphys_down = unphysical_points-1;
        unphys_down = unphys_down(unphys_down>=1);
        near_unphysical = unique([unphys_down; unphys_up]);
        plot_points = setdiff(near_unphysical, unphysical_points);

        not_unphysical = setdiff(1:length(phi_vals), unphysical_points);

        x = rad2deg(phi_vals);
        y = nan(size(E_n));
        y(plot_points) = E_n(plot_points);
        x = x(not_unphysical);
        y = y(not_unphysical);
        if use_newcolours
            plot(x, y, ':','Color',newcolors(counter), 'HandleVisibility','off');
        else
            plot(x, y, ':', 'HandleVisibility','off');
        end
    end

end

xlim([0,180])

if use_rel_E
    ylim([-1.1,1.1]);
end

% newcolours_cell_array_rgb = hex2rgb(cellstr(erase(newcolors(1:length(analysis_output.phi_at_min)),'#')))/255;
newcolours_cell_array_rgb = hex2rgb(cellstr(newcolors(1:length(analysis_output.phi_at_min))));

if plot_local_minima
%     scatter(rad2deg(analysis_output.phi_at_min), ...
%         (analysis_output.E_all_min(1,:)-E_unwrapped(:,1)')./scaling,50,...
%         'p','MarkerEdgeColor',p1.Color, 'MarkerFaceColor',p1.Color,...
%         'HandleVisibility','off')
    scatter(rad2deg(analysis_output.phi_loc_min), ...
        (analysis_output.E_loc_min(1,:)-E_unwrapped(:,1)')./scaling,50,...
        newcolours_cell_array_rgb,'filled','o',...
        'HandleVisibility','off')
    % also put dots at the endpoints if they're minima
    counter = 0;
    for jj = slice
        counter = counter+1;
        if analysis_output.zero_phi_is_min(counter)
            scatter(rad2deg(min(phi_vals)), ...
                (analysis_output.E_low_phi(counter)-E_unwrapped(counter,1)')./scaling(counter),50,...
                newcolours_cell_array_rgb(counter,:),'filled','o',...
                'HandleVisibility','off')
        end
        if analysis_output.pi_phi_is_min(counter)&&~dont_plot_max(counter)
            scatter(rad2deg(max(phi_vals)), ...
                (analysis_output.E_high_phi(counter)-E_unwrapped(counter,1)')./scaling(counter),50,...
                newcolours_cell_array_rgb(counter,:),'filled','o',...
                'HandleVisibility','off')
        end
    end
end

if plot_minima
%     scatter(rad2deg(analysis_output.phi_at_min), ...
%         (analysis_output.E_all_min(1,:)-E_unwrapped(:,1)')./scaling,50,...
%         'p','MarkerEdgeColor',p1.Color, 'MarkerFaceColor',p1.Color,...
%         'HandleVisibility','off')
    scatter(rad2deg(analysis_output.phi_at_min), ...
        (analysis_output.E_all_min(1,:)-E_unwrapped(:,1)')./scaling,300,...
        newcolours_cell_array_rgb,...
        'filled','p',...
        'HandleVisibility','off')
end

if include_legend
    legend
end

end