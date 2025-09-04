function [x1,y1] = plot_cases_phi(ps, outputs_analyse, type, varargin)

args = varargin;
nargs = numel(args);
plot_specific_points = false;
plot_local_minima = false;
markertype = '-';
markertype_local_min = ':';
use_newcolours = false;
handle_visibility = 'off';
UseLocalMinName = false;
plot_stars = false;
phi_limit = pi;
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
    elseif strcmpi(args{k},'markertype')
        markertype = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'markertype_local_min')||strcmpi(args{k},'MarkertypeLocalMin')
        markertype_local_min = args{k+1};
        % we want to only use this if we're actually plotting local minima
        assert(plot_local_minima);
        k = k+1;
    elseif strcmpi(args{k},'PlotLocalMinima')||strcmpi(args{k},'Plot_Local_Minima')
        plot_local_minima = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'Colour')||strcmpi(args{k},'Color')
        colour = args{k+1};
        use_newcolours = true;
        k = k+1;
    elseif strcmpi(args{k},'PlotStarsLimit')||strcmpi(args{k},'Plot_Stars_Limit')
        phi_limit = args{k+1};
        plot_stars = true;
        k = k+1;
    elseif strcmpi(args{k},'UseLocalMinName')
        UseLocalMinName = args{k+1};
        if UseLocalMinName
            handle_visibility = 'on';
        end
%         use_newcolours = false;
        k = k+1;
    else
        error('intput not recognised!');
    end
    k = k+1;
end


switch type
    case 4
        x1 = ps.kappa_slice./ps.kD_slice;
        displayname = sprintf('$\\epsilon = %0.3g$', ps.epsilon_slice(1));
    case 1
        x1 = -ps.kappa_slice./ps.epsilon_slice;
        displayname = sprintf('$k_D = %0.3g$', ps.kD_slice(1));
    case 2
        x1 = -ps.kappa_slice./ps.epsilon_slice;
        displayname = sprintf('$k_D = %0.3g$', ps.kD_slice(1));
    case 3
        x1 = ps.kappa_slice./ps.kD_slice;
        displayname = sprintf('$\\epsilon = %0.3g$', ps.epsilon_slice(1));
    case 5
        x1 = -ps.kD_slice./ps.epsilon_slice;
        displayname = sprintf('$\\kappa = %0.3g$', ps.kappa_slice(1));
    case 6
        x1 = -ps.kD_slice./ps.epsilon_slice;
        displayname =  sprintf('$\\kappa = %0.3g$', ps.kappa_slice(1));
    case 7
        x1 = ps.sigma_slice;
        displayname = 'true global minimum';
    case 8
        x1 = ps.zeta_slice;
        displayname = 'true global minimum';
end

y1 = outputs_analyse.phi_at_min;
[x,y] = sort_by_x(x1,y1);
if use_newcolours
    v1 = plot(x, y,markertype,'displayname', displayname, 'Color',colour);
else
    v1 = plot(x, y,markertype,'displayname', displayname);
end

if plot_local_minima
    if isfield(outputs_analyse, 'phi_local_min')
        y2 = nan(size(x));
        for ii=1:length(outputs_analyse.phi_local_min)
            local_mins = outputs_analyse.phi_local_min{ii};
            if ~isempty(local_mins)
                y2(ii) = local_mins(1);
            else
                y2(ii) = nan;
            end
        end
    elseif isfield(outputs_analyse, 'phi_loc_min')
        y2 = outputs_analyse.phi_loc_min;
%         [x,y2] = sort_by_x(x,y2);
    end
    [x,y2] = sort_by_x(x1,y2);
    if use_newcolours
        v2 = plot(x, y2,markertype_local_min,'HandleVisibility',handle_visibility,...
            'Color',colour,'DisplayName', 'local minimum for small phi');
    else
        v2 = plot(x, y2,markertype_local_min,'HandleVisibility',handle_visibility, ...
            'Displayname','local minimum for small phi', 'Color', "#FE6100");
    end



end

if plot_stars
    % find the point where phi gets to a maximum, and plot a star at that
    % point and one with smaller phi next to it
    index = find(y<phi_limit,1)-1;
    if ~(index<1||index>length(x))
        [val_smaller, index_smaller] = min(y(index-1:2:index+1));
        if index_smaller==1
            index_smaller = index-1;
        else
            index_smaller = index+1;
        end
        % now plot stars at those indices
        scatter(x([index, index_smaller]), y([index, index_smaller]),300,...
            'filled','p','HandleVisibility','off', 'MarkerEdgeColor',colour, ...
            'MarkerFaceColor',colour)
    end
end


end