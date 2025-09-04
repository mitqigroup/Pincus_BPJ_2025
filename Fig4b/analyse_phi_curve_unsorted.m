function output = analyse_phi_curve_unsorted(E_all,data_structure,slice,phi_vals,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% can also try smoothdata, islocalmax, islocalmin (auto matlab functions)

% plot_curves = options.plot_curves;

% parse inputs and set defaults
args = varargin{1};
nargs = numel(args);
plot_curves = false;
ignore_near_nan = false;
remove_unphysical = false;
clean_data = false;
% check_fake_maxima = false;
k = 1;
while k<=nargs
    % if we get a figure handle, use it
    if strcmpi(args{k},'PlotCurves')||strcmpi(args{k},'Plot_Curves')
        plot_curves = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'xlabel')
        xaxis_name = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'IgnoreNearNaN')
        ignore_near_nan = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'RemoveUnphysical')
        remove_unphysical = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'CleanData')||strcmpi(args{k},'Clean_Data')
        clean_data = args{k+1};
        k = k+1;
    else
        error('intput not recognised!');
    end
    k = k+1;
end

%% analyse data from each, we want local minima and maxima and end points.
% for jj = 1:length(parameter_set)
% plot_curves = true;
% plot_curves = false;
[~, lowest_phi_index] = min(phi_vals);
[~, highest_phi_index] = max(phi_vals);
phi_loc_min = nan(size(slice));
E_loc_min = nan(size(slice));
phi_loc_max = nan(size(slice));
E_loc_max = nan(size(slice));
for counter = 1:length(slice)
    jj = slice(counter);
    local_min = false;
    local_max = false;

    E_all_this = squeeze(E_all(:,:,jj));
    Sigma_phi = squeeze(data_structure.Sigma_vals_nonlinear(:,jj));
    aB_phi = squeeze(data_structure.alpha_B_vals_nonlinear(:,jj));

    unphysical_removed = [];
    if remove_unphysical
        % remove all the unphysical points. First remove all of those with
        % a clearly unphysical phi_end, and store the locations
        if isfield(data_structure, 'phi_end_nonlinear')
            phi_end_this = data_structure.phi_end_nonlinear(:,jj);
            unphysical_removed = cat(1,unphysical_removed, find(phi_end_this>0.1));
        end
        % now we want to also remove ones where there's an unusually large
        % jump in alpha_B, since this says that it's not on the same branch
        % as the other ones
        diff_aB = diff(aB_phi);
        % check points where we're 5 times the mean, and large phi
        new_remove = find((abs(diff_aB'/mean(diff_aB))>5).*(phi_vals(2:end)>2.6116));
%         phi_vals(2:end)>2.6116
        % also remove all points between the two limits
        new_remove = min(new_remove):max(new_remove);
        unphysical_removed = unique(cat(1,unphysical_removed, new_remove'));
        E_all_this(:,unphysical_removed) = nan;
    end

    if clean_data
        % this is a general cleaning and smoothing of the data. We
        % eliminate local C0 discontinuities in both alpha_B and the total
        % energy, then excise those points and store them for later
        % plotting. 

    end

%     E = squeeze(E_all(1,:,jj));
%     E_ad = squeeze(E_all(2,:,jj));
%     E_bend_A = squeeze(E_all(5,:,jj));
%     E_bend_B = squeeze(E_all(6,:,jj));
%     E_stretch_A = squeeze(E_all(3,:,jj));
%     E_stretch_B = squeeze(E_all(4,:,jj));
    E = E_all_this(1,:);
    E_ad = E_all_this(2,:);
    E_bend_A = E_all_this(5,:);
    E_bend_B = E_all_this(6,:);
    E_stretch_A = E_all_this(3,:);
    E_stretch_B = E_all_this(4,:);
    diff_E = diff(E);

    if plot_curves
        figure();
        hold on
        plot(phi_vals,E)
    end

    % find points where diff changes sign, these are stationary points.
    locs_pos = [];
    locs_neg = [];
    for ii = 2:length(diff_E)
        if ~(ignore_near_nan&&any(isnan(E(max([ii-2,1]):min([ii+2,length(E)])))))
            if diff_E(ii)<0&&diff_E(ii-1)>0
                locs_pos = [locs_pos, ii];
            elseif diff_E(ii)>0&&diff_E(ii-1)<0
                locs_neg = [locs_neg, ii];
            end
        end
    end
    
    if plot_curves
        scatter(phi_vals(locs_neg), E(locs_neg), 'pk')
        scatter(phi_vals(locs_pos), E(locs_pos), 'sk')
    end

    % there should be at most one local maxima and minima, if it exists get
    % the exact location via a quadratic fit
    if ~isempty(locs_pos)&&length(locs_pos)==1
        local_max = true;
        max_range =(locs_pos(1)-1):(locs_pos(1)+1);
        phi_near_max = phi_vals(max_range);
        E_near_max = E(max_range);
        [E_loc_max(counter), phi_loc_max(counter), polycoeff] = get_polyfit(phi_near_max, E_near_max);
        if plot_curves
            x = linspace(phi_near_max(1),phi_near_max(end));
            plot(x, polyval(polycoeff, x))
        end
    elseif ~isempty(locs_pos)
        fprintf('Multiple local maxima for slice=%i, investigate \n',jj);
        local_max = true;
    end
    if ~isempty(locs_neg)&&length(locs_neg)==1
        local_min = true;
        min_range =(locs_neg(1)-1):(locs_neg(1)+1);
        phi_near_min = phi_vals(min_range);
        E_near_min = E(min_range);
        [E_loc_min(counter), phi_loc_min(counter), polycoeff] = get_polyfit(phi_near_min, E_near_min);
        E_ad_lm = polyval(polyfit(phi_near_min, E_ad(min_range), 2), phi_loc_min(counter));
        E_bendA_lm = polyval(polyfit(phi_near_min, E_bend_A(min_range), 2), phi_loc_min(counter));
        E_bendB_lm = polyval(polyfit(phi_near_min, E_bend_B(min_range), 2), phi_loc_min(counter));
        E_stretchA_lm = polyval(polyfit(phi_near_min, E_stretch_A(min_range), 2), phi_loc_min(counter));
        E_stretchB_lm = polyval(polyfit(phi_near_min, E_stretch_B(min_range), 2), phi_loc_min(counter));
        Sigma_lm = polyval(polyfit(phi_near_min, Sigma_phi(min_range), 2), phi_loc_min(counter));
        if plot_curves
            x = linspace(phi_near_min(1),phi_near_min(end));
            plot(x, polyval(polycoeff, x))
        end
        loc_min_index = locs_neg(1);
    elseif ~isempty(locs_neg)
        fprintf('Multiple local minima for slice=%i, investigate \n',jj);
        local_min = true;
        % get the local minima with the smallest phi
        [~,loc_min_index] = min(phi_vals(locs_neg));
        min_range =(locs_neg(loc_min_index)-1):(locs_neg(loc_min_index)+1);
        phi_near_min = phi_vals(min_range);
        E_near_min = E(min_range);
        [E_loc_min(counter), phi_loc_min(counter), polycoeff] = get_polyfit(phi_near_min, E_near_min);
        E_ad_lm = polyval(polyfit(phi_near_min, E_ad(min_range), 2), phi_loc_min(counter));
        E_bendA_lm = polyval(polyfit(phi_near_min, E_bend_A(min_range), 2), phi_loc_min(counter));
        E_bendB_lm = polyval(polyfit(phi_near_min, E_bend_B(min_range), 2), phi_loc_min(counter));
        E_stretchA_lm = polyval(polyfit(phi_near_min, E_stretch_A(min_range), 2), phi_loc_min(counter));
        E_stretchB_lm = polyval(polyfit(phi_near_min, E_stretch_B(min_range), 2), phi_loc_min(counter));
        Sigma_lm = polyval(polyfit(phi_near_min, Sigma_phi(min_range), 2), phi_loc_min(counter));
        if plot_curves
            x = linspace(phi_near_min(1),phi_near_min(end));
            plot(x, polyval(polycoeff, x))
        end
    end

    E_low_phi(counter) = E(lowest_phi_index);
    E_high_phi(counter) = E(highest_phi_index);
    [min_val,min_index] = min(E);
    [max_val,max_index] = max(E);

    % check if the endpoints are minima, different branches for ascending
    % or descending phi
    % first the lowest phi
    diff_small = diff(squeeze(E_all(1,1:2,jj)));
    diff_large = diff(squeeze(E_all(1,end-1:end,jj)));
    if lowest_phi_index==1
        if diff_small>0
            zero_phi_is_min(counter) = true;
        else
            zero_phi_is_min(counter) = false;
        end
        if diff_large<0
            pi_phi_is_min(counter) = true;
        else
            pi_phi_is_min(counter) = false;
        end
    else
        if diff_large>0
            zero_phi_is_min(counter) = false;
        else
            zero_phi_is_min(counter) = true;
        end
        if diff_small<0
            pi_phi_is_min(counter) = false;
        else
            pi_phi_is_min(counter) = true;
        end
    end

    % go through the possibilities
    % first is no local minima/maxima, min is at high phi, that implies
    % fully attached (curve type 1)
    if all(isnan(E))
        curve_type(counter) = 0;
        phi_at_min(counter) = nan;
%         E_all_min(:,counter) = squeeze(E_all(:,highest_phi_index,jj));
        E_all_min(:,counter) = nan;
        phi_at_max(counter) = nan;
        Sigma_at_min(counter) = nan;
        fprintf('All nans lol \n',jj);
    elseif ~local_max&&~local_min&&min_index==highest_phi_index
        curve_type(counter) = 1;
        phi_at_min(counter) = phi_vals(highest_phi_index);
%         E_all_min(:,counter) = squeeze(E_all(:,highest_phi_index,jj));
        E_all_min(:,counter) = E_all_this(:,highest_phi_index);
        phi_at_max(counter) = phi_vals(lowest_phi_index);
        Sigma_at_min(counter) = Sigma_phi(highest_phi_index);
    % second is local minima and maxima, but min is at highest phi
    elseif local_max&&local_min&&min_index==highest_phi_index
        curve_type(counter) = 2;
        phi_at_min(counter) = phi_vals(highest_phi_index);
%         E_all_min(:,counter) = squeeze(E_all(:,highest_phi_index,jj));
        E_all_min(:,counter) = E_all_this(:,highest_phi_index);
        phi_at_max(counter) = phi_vals(lowest_phi_index);
        Sigma_at_min(counter) = Sigma_phi(highest_phi_index);
    % third is local minima and maxima, but min is at local min
    elseif local_max&&local_min&&min_index==loc_min_index
        curve_type(counter) = 3;
        phi_at_min(counter) = phi_loc_min(counter);
        E_all_min(:,counter) = [E_loc_min(counter),...
            E_ad_lm,E_stretchA_lm, E_stretchB_lm,E_bendA_lm,E_bendB_lm];
        phi_at_max(counter) = phi_vals(max_index);
        Sigma_at_min(counter) = Sigma_lm;
    % fourth is local minima but no local maxima, min is at local min (might not
    % exist)
    elseif ~local_max&&local_min&&min_index==loc_min_index
        curve_type(counter) = 7;
        phi_at_min(counter) = phi_loc_min(counter);
        E_all_min(:,counter) = [E_loc_min(counter),...
            E_ad_lm,E_stretchA_lm, E_stretchB_lm,E_bendA_lm,E_bendB_lm];
        phi_at_max(counter) = phi_vals(max_index);
        Sigma_at_min(counter) = Sigma_lm;
    % fifth is no local minima but local maxima, min is at 0
    elseif local_max&&~local_min&&min_index==lowest_phi_index
        curve_type(counter) = 6;
        phi_at_min(counter) = phi_vals(lowest_phi_index);
%         E_all_min(:,counter) = squeeze(E_all(:,lowest_phi_index,jj));
        E_all_min(:,counter) = E_all_this(:,lowest_phi_index);
        phi_at_max(counter) = phi_vals(max_index);
        Sigma_at_min(counter) = Sigma_phi(lowest_phi_index);
    % sixth is no local minima or maxima, min is at 0
    elseif ~local_max&&~local_min&&min_index==lowest_phi_index
        curve_type(counter) = 5;
        phi_at_min(counter) = phi_vals(lowest_phi_index);
%         E_all_min(:,counter) = squeeze(E_all(:,lowest_phi_index,jj));
        E_all_min(:,counter) = E_all_this(:,lowest_phi_index);
        phi_at_max(counter) = phi_vals(highest_phi_index);
        Sigma_at_min(counter) = Sigma_phi(lowest_phi_index);
    % seventh is unwrapped, but local minima exists. We don't care about
    % maxima for this one. Likely degenerate (nope!)
    elseif local_min&&min_index==lowest_phi_index
        curve_type(counter) = 4;
        phi_at_min(counter) = phi_vals(lowest_phi_index);
%         E_all_min(:,counter) = squeeze(E_all(:,lowest_phi_index,jj));
        E_all_min(:,counter) = E_all_this(:,lowest_phi_index);
        phi_at_max(counter) = phi_vals(highest_phi_index);
        Sigma_at_min(counter) = Sigma_phi(lowest_phi_index);
    % this one is almost certainly degenerate
    elseif local_max&&min_index==lowest_phi_index
        curve_type(counter) = 8;
        phi_at_min(counter) = phi_vals(lowest_phi_index);
%         E_all_min(:,counter) = squeeze(E_all(:,lowest_phi_index,jj));
        E_all_min(:,counter) = E_all_this(:,lowest_phi_index);
        phi_at_max(counter) = phi_vals(highest_phi_index);
        Sigma_at_min(counter) = Sigma_phi(lowest_phi_index);
    else
        curve_type(counter) = 0;
        phi_at_min(counter) = phi_vals(min_index);
%         E_all_min(:,counter) = squeeze(E_all(:,min_index,jj));
        E_all_min(:,counter) = E_all_this(:,min_index);
        phi_at_max(counter) = phi_vals(max_index);
        Sigma_at_min(counter) = Sigma_phi(min_index);
        fprintf('You have a new curve type for slice=%i, investigate!! \n',jj);
    end

    % get local maxima and minima locations, put in cell array
    local_minima_locations{counter} = locs_neg;
    local_maxima_locations{counter} = locs_pos;
    unphysical_points{counter} = unphysical_removed;

end

output.phi_at_min = phi_at_min;
output.phi_at_max = phi_at_max;
output.curve_type = curve_type;
output.E_all_min = E_all_min;
output.phi_loc_min = phi_loc_min;
output.E_loc_min = E_loc_min;
output.Sigma_at_min = Sigma_at_min;
output.local_minima_locations = local_minima_locations;
output.local_maxima_locations = local_maxima_locations;
output.unphysical_points = unphysical_points;
output.zero_phi_is_min = zero_phi_is_min;
output.pi_phi_is_min = pi_phi_is_min;
output.E_low_phi = E_low_phi;
output.E_high_phi = E_high_phi;

end





