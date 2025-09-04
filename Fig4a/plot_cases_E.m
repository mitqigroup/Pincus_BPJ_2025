function [x,y] = plot_cases_E(ps, outputs_analyse, type, varargin)

% parse inputs and set defaults
args = varargin;
nargs = numel(args);
plot_specific_points = false;
k = 1;
while k<=nargs
    % if we get a figure handle, use it
    if strcmpi(args{k},'PlotSpecificPoints')||strcmpi(args{k},'Plot_Specific_Points')
        plotting_points = args{k+1};
        plot_specific_points = true;
        % this should be a cell array of specific points to plot, with the
        % same size as the slices
        assert(length(slices)==length(plotting_points), 'Slices not same as plotting points');
        k = k+1;
    else
        error('intput not recognised!');
    end
    k = k+1;
end

% markertype = '.-';
markertype = '-';

if isfield(outputs_analyse,'all_global_minima')
    s = cell2mat(outputs_analyse.all_global_minima);
    switch type
        case 1
            x = -ps.kappa_slice./ps.epsilon_slice;
    
            c = {s.E_bendA_min};
            ybA = cell2mat(c);
            c = {s.E_bendB_min};
            ybB = cell2mat(c);
            c = {s.E_adhesion_min};
            yad = cell2mat(c);
            y = -(ybA+ybB)./yad;
            y(outputs_analyse.phi_at_min<0.01) = nan;
    
            displayname = sprintf('$k_D = %0.3g$', ps.kD_slice(1));
        case 2
            x = -ps.kappa_slice./ps.epsilon_slice;
    
            c = {s.E_bendA_min};
            ybA = cell2mat(c);
            c = {s.E_bendB_min};
            ybB = cell2mat(c);
            c = {s.E_adhesion_min};
            yad = cell2mat(c);
            y = -(ybA+ybB)./yad;
            y(outputs_analyse.phi_at_min<0.01) = nan;
            
            displayname = sprintf('$k_D = %0.3g$', ps.kD_slice(1));
        case 3
            E_unwrapped = ps.kD_slice./2.*ps.alpha_i_slice.^2.*ps.d_slice.^2./(1+ps.alpha_i_slice);
            x = ps.kappa_slice./ps.kD_slice;
    
            c = {s.E_bendA_min};
            ybA = cell2mat(c);
            c = {s.E_bendB_min};
            ybB = cell2mat(c);
            c = {s.E_stretchA_min};
            ysA = cell2mat(c);
            c = {s.E_stretchB_min};
            ysB = cell2mat(c);
            y = -(ybA+ybB)./(ysA+ysB-E_unwrapped);
            y(outputs_analyse.phi_at_min<0.01) = nan;

            displayname = sprintf('$\\epsilon = %0.3g$', ps.epsilon_slice(1));
        case 4
            E_unwrapped = ps.kD_slice./2.*ps.alpha_i_slice.^2.*ps.d_slice.^2./(1+ps.alpha_i_slice);
            x = ps.kappa_slice./ps.kD_slice;
    
            c = {s.E_bendA_min};
            ybA = cell2mat(c);
            c = {s.E_bendB_min};
            ybB = cell2mat(c);
            c = {s.E_stretchA_min};
            ysA = cell2mat(c);
            c = {s.E_stretchB_min};
            ysB = cell2mat(c);
            y = -(ybA+ybB)./(ysA+ysB-E_unwrapped);
            y(outputs_analyse.phi_at_min<0.01) = nan;

            displayname = sprintf('$\\epsilon = %0.3g$', ps.epsilon_slice(1));
        case 5
            E_unwrapped = ps.kD_slice./2.*ps.alpha_i_slice.^2.*ps.d_slice.^2./(1+ps.alpha_i_slice);
            x = -ps.kD_slice./ps.epsilon_slice;
    
            c = {s.E_stretchA_min};
            ysA = cell2mat(c);
            c = {s.E_stretchB_min};
            ysB = cell2mat(c);
            c = {s.E_adhesion_min};
            yad = cell2mat(c);
            y = -(ysA+ysB-E_unwrapped)./(yad);
            y(outputs_analyse.phi_at_min<0.01) = nan;

            displayname = sprintf('$\\kappa = %0.3g$', ps.kappa_slice(1));
        case 6
            E_unwrapped = ps.kD_slice./2.*ps.alpha_i_slice.^2.*ps.d_slice.^2./(1+ps.alpha_i_slice);
            x = -ps.kD_slice./ps.epsilon_slice;
    
            c = {s.E_stretchA_min};
            ysA = cell2mat(c);
            c = {s.E_stretchB_min};
            ysB = cell2mat(c);
            c = {s.E_adhesion_min};
            yad = cell2mat(c);
            y = -(ysA+ysB-E_unwrapped)./(yad);
            y(outputs_analyse.phi_at_min<0.01) = nan;

            displayname = sprintf('$\\kappa = %0.3g$', ps.kappa_slice(1));
    end
    [x,y] = sort_by_x(x,y);
    v1 = plot(x, y,markertype,'displayname',displayname);
else
    % old version
    switch type
        case 4
            E_unwrapped = ps.kD_slice./2.*ps.alpha_i_slice.^2.*ps.d_slice.^2./(1+ps.alpha_i_slice);
            x = ps.kappa_slice./ps.kD_slice;
            y = (sum(outputs_analyse.E_all_min(5:6,:),1))./(sum(outputs_analyse.E_all_min(3:4,:),1)-E_unwrapped');
            y(outputs_analyse.phi_at_min<0.01) = nan;
            [x,y] = sort_by_x(x,y);
            v1 = plot(x, y,markertype,'displayname', sprintf('$\\epsilon = %0.3g$', ps.epsilon_slice(1)));
        case 1
            x = -ps.kappa_slice./ps.epsilon_slice;
            y = -(sum(outputs_analyse.E_all_min(5:6,:),1))./(outputs_analyse.E_all_min(2,:));
            y(outputs_analyse.phi_at_min<0.01) = nan;
            [x,y] = sort_by_x(x,y);
            v1 = plot(x, y,markertype,'displayname', sprintf('$k_D = %0.3g$', ps.kD_slice(1)));
        case 2
            x = -ps.kappa_slice./ps.epsilon_slice;
            y = -(sum(outputs_analyse.E_all_min(5:6,:),1))./(outputs_analyse.E_all_min(2,:));
            y(outputs_analyse.phi_at_min<0.01) = nan;
            [x,y] = sort_by_x(x,y);
            v1 = plot(x, y,markertype,'displayname', sprintf('$k_D = %0.3g$', ps.kD_slice(1)));
        case 3
            E_unwrapped = ps.kD_slice./2.*ps.alpha_i_slice.^2.*ps.d_slice.^2./(1+ps.alpha_i_slice);
            x = ps.kappa_slice./ps.kD_slice;
            y = (sum(outputs_analyse.E_all_min(5:6,:),1))./(sum(outputs_analyse.E_all_min(3:4,:),1)-E_unwrapped');
            y(outputs_analyse.phi_at_min<0.01) = nan;
            [x,y] = sort_by_x(x,y);
            v1 = plot(x, y,markertype,'displayname',sprintf('$\\epsilon = %0.3g$', ps.epsilon_slice(1)));
        case 5
            E_unwrapped = ps.kD_slice./2.*ps.alpha_i_slice.^2.*ps.d_slice.^2./(1+ps.alpha_i_slice);
            x = -ps.kD_slice./ps.epsilon_slice;
            y = -(sum(outputs_analyse.E_all_min(3:4,:),1)-E_unwrapped')./(outputs_analyse.E_all_min(2,:));
            y(outputs_analyse.phi_at_min<0.01) = nan;
            [x,y] = sort_by_x(x,y);
            v1 = plot(x, y,markertype,'displayname',sprintf('$\\kappa = %0.3g$', ps.kappa_slice(1)));
        case 6
            E_unwrapped = ps.kD_slice./2.*ps.alpha_i_slice.^2.*ps.d_slice.^2./(1+ps.alpha_i_slice);
            x = -ps.kD_slice./ps.epsilon_slice;
            y = -(sum(outputs_analyse.E_all_min(3:4,:),1)-E_unwrapped')./(outputs_analyse.E_all_min(2,:));
            y(outputs_analyse.phi_at_min<0.01) = nan;
            [x,y] = sort_by_x(x,y);
            v1 = plot(x, y,markertype,'displayname',sprintf('$\\kappa = %0.3g$', ps.kappa_slice(1)));
    end
end

end