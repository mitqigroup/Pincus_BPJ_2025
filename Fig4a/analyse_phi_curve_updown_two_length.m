function output = analyse_phi_curve_updown_two_length(ds_middle, ds_down,slice,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

warning('off', 'SPLINES:CHCKXYWP:NaNs');

% parse inputs and set defaults
args = varargin;
nargs = numel(args);
plot_curves = false;
movmedian_range = 3;
movmedian_threshold = 0.01;
weights_range = 3;
weights_val = 1e3;
use_weights = false;
N_spline_points = 10000;
k = 1;
while k<=nargs
    % if we get a figure handle, use it
    if strcmpi(args{k},'PlotCurves')||strcmpi(args{k},'Plot_Curves')
        plot_curves = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'MovmedianRange')||strcmpi(args{k},'Movmedian_Range')
        movmedian_range = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'MovmedianThreshold')||strcmpi(args{k},'Movmedian_Threshold')
        movmedian_threshold = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'WeightsRange')||strcmpi(args{k},'Weights_Range')
        weights_range = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'UseWeights')||strcmpi(args{k},'Use_Weights')
        use_weights = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'WeightsEndValues')||strcmpi(args{k},'Weights_End_Values')
        weights_val = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'NSplinePoints')||strcmpi(args{k},'NumberSplinePoints')
        N_spline_points = args{k+1};
        k = k+1;
    else
        error('intput not recognised!');
    end
    k = k+1;
end

if ~ds_middle.is_sorted
    error('for this script, data must be sorted');
end

%% analyse data from each, we want local minima and maxima and end points.
for counter = 1:length(slice)
    jj = slice(counter);

    E_middle = ds_middle.E(jj,:);
    E_down = ds_down.E(jj,:);
    remove_points = false(size(E_middle));

    E_std = std([E_middle;E_down])./mean([E_middle;E_down]);
    E = E_middle;
    E_removed = E;
    remove_points(abs(E_std)>1e-3) = true;
    E_removed(remove_points) = nan;
    [outliers,L,U,C] = isoutlier(E_removed, 'movmedian',movmedian_range, ...
        'ThresholdFactor',movmedian_threshold);
    outliers(1:floor(movmedian_range/2)) = false;
    outliers(end-floor(movmedian_range/2):end) = false;
    % sometimes the final point is very incorrect, check slope just in case
    diff_E = diff(E);
    if abs(diff_E(end))>mean(abs(diff_E))
        remove_points(end) = true;
    end
    remove_points(outliers) = true;
    
    xx = linspace(0,pi,N_spline_points);
%     fprintf('here %d \n', counter);
    E_spline = fit_spline(ds_middle.phi_vals, E, xx, remove_points, ...
        'UseWeights', use_weights,'WeightsRange', weights_range, ...
        'WeightsEndValues', weights_val);
    % this should be extended to fit all the possible parameters, alpha,
    % Sigma, E_bend etc.
    
    [local_minima, min_promenance] = islocalmin(E_spline);
    [local_maxima, max_promenance] = islocalmax(E_spline);
    % check if first and last points are global minima
    phi_zero_is_minima = false;
    phi_pi_is_minima = false;
    if E_spline(1)<E_spline(2)
        phi_zero_is_minima  = true;
    end
    if E_spline(end)<E_spline(end-1)
        phi_pi_is_minima = true;
    end
    [global_min, global_min_index] = min(E_spline);
    [global_max, global_max_index] = max(E_spline);

    if plot_curves
        figure();
        hold on
        plot(ds_middle.phi_vals, E_middle);
        plot(ds_middle.phi_vals(~remove_points), E_middle(~remove_points));
        plot(xx, E_spline);
        plot(xx(local_minima), E_spline(local_minima), 'ro');
        plot(xx(local_maxima), E_spline(local_maxima), 'rs');
        plot(ds_middle.phi_vals, L, 'k:', 'LineWidth',0.5);
        plot(ds_middle.phi_vals, U, 'k:', 'LineWidth',0.5);
%         yyaxis right
%         plot(ds_middle.phi_vals,psi_end_middle);
%         plot(ds_middle.phi_vals,psi_end_down);
%         yyaxis left

    end

    % get all minima
    yy = fit_spline(ds_middle.phi_vals, ds_middle.E(jj,:), xx, remove_points);
    E_min_global(counter) = yy(global_min_index);
    yy = fit_spline(ds_middle.phi_vals, ds_middle.aA(jj,:), xx, remove_points);
    aA_min_global(counter) = yy(global_min_index);
    yy = fit_spline(ds_middle.phi_vals, ds_middle.aB(jj,:), xx, remove_points);
    aB_min_global(counter) = yy(global_min_index);
    phi_local_min{counter} = xx(local_minima);
    E_local_min{counter} = E_spline(local_minima);
    phi_local_max{counter} = xx(local_maxima);
    E_local_max{counter} = E_spline(local_maxima);


    % classify curves
    if all(isnan(E))
        curve_type(counter) = 0;
        fprintf('All nans lol \n',jj);
    elseif all(~local_minima)&&global_min_index==length(E_spline)
        curve_type(counter) = 1;
    elseif any(local_minima)&&global_min_index==length(E_spline)
        curve_type(counter) = 2;
    elseif any(local_minima)&&global_min_index==find(local_minima,1)&&phi_pi_is_minima
        curve_type(counter) = 3;
    elseif any(local_minima)&&global_min_index==1&&(E(end)>E(end-1))
        curve_type(counter) = 4;
    elseif all(~local_minima)&&global_min_index==1&&global_max_index==length(E_spline)
        curve_type(counter) = 5;
    elseif all(~local_minima)&&global_min_index==1&&global_max_index==find(local_maxima,1)
        curve_type(counter) = 6;
    elseif any(local_minima)&&global_min_index==find(local_minima,1)&&~phi_pi_is_minima&&length(find(local_minima,1))==1
        curve_type(counter) = 3;
    elseif global_min_index==1&&(E(end)<E(end-1))
        curve_type(counter) = 6;
    else
        curve_type(counter) = 0;
        fprintf('You have a new curve type for slice=%i, investigate!! \n',jj);
    end

    % get local maxima and minima locations, put in cell array
%     local_minima_locations{counter} = local_minima;
%     local_maxima_locations{counter} = local_maxima;
    removed_points{counter} = remove_points;
    zero_phi_is_min(counter) = phi_zero_is_minima;
    pi_phi_is_min(counter) = phi_pi_is_minima;
    phi_at_min(counter) = xx(global_min_index);
    phi_at_max(counter) = xx(global_max_index);

end

output.phi_at_min = phi_at_min;
output.phi_at_max = phi_at_max;
output.curve_type = curve_type;
output.unphysical_points = removed_points;
output.zero_phi_is_min = zero_phi_is_min;
output.pi_phi_is_min = pi_phi_is_min;

output.phi_local_min = phi_local_min;
output.E_local_min = E_local_min;
output.phi_local_max = phi_local_max;
output.E_local_max = E_local_max;

output.E_min_global = E_min_global;
output.aA_min_global = aA_min_global;
output.aB_min_global = aB_min_global;

end





