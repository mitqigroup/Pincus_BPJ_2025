function fig_handle = plot_minima_xy(x,y, varargin)

% parse inputs and set defaults
args = varargin;
nargs = numel(args);
k = 1;
fig_input = false;
use_newcolours = false;
update_yticks = false;
update_ytick_labels = false;
update_ylimits = false;
xaxis_name = '$\phi$';
yaxis_name = '$\Delta E$';
displayname = '$y$';
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
    elseif strcmpi(args{k},'ColourOrder')||strcmpi(args{k},'ColorOrder')
        newcolors = args{k+1};
        use_newcolours = true;
        k = k+1;
    elseif strcmpi(args{k},'yTicks')
        ytick_locs = args{k+1};
        update_yticks = true;
        k = k+1;
    elseif strcmpi(args{k},'yTickLabels')
        ytick_labels = args{k+1};
        update_ytick_labels = true;
        k = k+1;
    elseif strcmpi(args{k},'ylim')||strcmpi(args{k},'ylimits')||strcmpi(args{k},'y_limits')
        y_limits = args{k+1};
        update_ylimits = true;
        k = k+1;
    elseif strcmpi(args{k},'DisplayName')||strcmpi(args{k},'Display_name')
        displayname = args{k+1};
        k = k+1;
    end
    k = k+1;
end

if fig_input
    fig_handle = figure(fig_handle_in,'Position',[400,100,800,600]);
else
    fig_handle = figure('Position',[400,100,800,600]);
end
axes1 = gca;
hold on
xlabel(xaxis_name)
ylabel(yaxis_name)
axes1.XScale = 'log';
if update_yticks
    yticks(ytick_locs)
%     yticks([0, pi/4, pi/2, 3*pi/4, pi])
end
if update_ytick_labels
    yticklabels(ytick_labels);
%     yticklabels({'$0$', '$\pi/4$', '$\pi/2$', '$3 \pi/4$', '$\pi$'})
end
if update_ylimits
    ylim(y_limits)
%     ylim([0,pi+0.01])
end
v1 = plot(x, y, 'displayname', displayname);
legend
if use_newcolours
    colororder(newcolors)
end
end