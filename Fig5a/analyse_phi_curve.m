function output = analyse_phi_curve(E_all,data_structure,slice,phi_vals,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% can also try smoothdata, islocalmax, islocalmin (auto matlab functions)

ds = data_structure;

% plot_curves = options.plot_curves;

% parse inputs and set defaults
args = varargin;
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

% if we don't have sorted data, use the old version of this script. The
% new version assumes that we have properly sorted data from now on, which
% is strictly ascending in phi
if ~data_structure.is_sorted
    output = analyse_phi_curve_unsorted(E_all,data_structure,slice,phi_vals,varargin);
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

    if plot_curves
        f1 = figure();
        hold on
        plot(phi_vals,E_all_this(1,:),'-')
    end

    unphysical_removed = [];
    if clean_data
        % x = rad2deg(ds.phi_vals');
        x = ds.phi_vals';
        % y = ds.alpha_B_vals_nonlinear(:,slice)-alpha_i_slice(1)-sigma_slice(1)*(1-cos(ds.phi_vals')).^2;
        y = ds.alpha_B_vals_nonlinear(:,jj)-ds.alpha_B_vals_nonlinear(1,jj);
        
        yf = nan(size(y));
        yb = nan(size(y));
        for ii=5:length(x)
            % get the last four points
            xi = x(ii-4:ii-1);
            yi = y(ii-4:ii-1);
            p = polyfit(xi,yi,2);
            yf(ii) = polyval(p, x(ii));
        end
        for ii=1:length(x)-4
            % get the last four points
            xi = x(ii+1:ii+4);
            yi = y(ii+1:ii+4);
            p = polyfit(xi,yi,2);
            yb(ii) = polyval(p, x(ii));
        end
        if plot_curves
            figure();
            hold on
            plot(x,y);
            plot(x,yf);
            plot(x,yb);
        end
        
        
        % ydf = abs(yf-y);
        % ydb = abs(yb-y);
        ydf = yf-y;
        ydb = yb-y;
%         ydf = ydf/mean(abs(ydf), 'omitnan');
%         ydb = ydb/mean(abs(ydb), 'omitnan');
        ydf = ydf/mean(abs(y), 'omitnan');
        ydb = ydb/mean(abs(y), 'omitnan');
        if plot_curves
            figure();
            hold on
            axes1 = gca;
            % axes1.YScale = 'log';
            plot(x, ydf);
            plot(x, ydb);
        end
        factor = 1e-2;
        i_ydf_pos = find(ydf>factor);
        i_ydb_pos = find(ydb>factor);
        i_ydf_neg = find(ydf<-factor);
        i_ydb_neg = find(ydb<-factor);
        
        if plot_curves
            figure();
            hold on
            plot(x,y);
            plot(x([i_ydb_pos;i_ydb_neg]), y([i_ydb_pos;i_ydb_neg]),'o');
            plot(x([i_ydf_pos;i_ydf_neg]), y([i_ydf_pos;i_ydf_neg]),'p');
        end
        
        % now for the forward one, there may be issues at high phi. Check if the
        % final point is wrong, then move backwards
        points_to_remove = [];
        check_point = length(ds.phi_vals);
        while any(i_ydf_pos==check_point)||any(i_ydf_neg==check_point)
            points_to_remove = [points_to_remove,check_point];
            check_point = check_point - 1;
        end
        % take out those points from i_ydf
        % points_to_remove
        i_ydf_pos = setdiff(i_ydf_pos,points_to_remove);
        i_ydf_neg = setdiff(i_ydf_neg,points_to_remove);
        
        % also take out all points from backwards which are larger than the new
        % largest forwards
        
        
        % get the union of both sets
        bf_unique = unique([i_ydb_pos; i_ydf_pos; i_ydb_neg; i_ydf_neg]);
        % get the smallest from the negative forwards
        min_f_neg = min(i_ydf_neg);
        % get the largest from the negative backwards
        max_b_neg = max(i_ydb_neg);
        % check that they're linked by the union. In other words, all values
        % between them are in the union. The difference should be empty!
        range_exclude = min_f_neg:max_b_neg;
        if isempty(setdiff(range_exclude, bf_unique))
            points_to_remove = [points_to_remove, range_exclude];
        end
        
        points_to_remove;
        
        if plot_curves
            figure();
            hold on
            E = ds.E_all_nonlinear(1,:,slice);
            plot(ds.phi_vals, E);
            E(points_to_remove) = nan;
            plot(ds.phi_vals, E, ':');
        end

        E_all_this(:,points_to_remove) = nan;

        unphysical_removed = unique(cat(1,unphysical_removed, points_to_remove'));

    end

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
%         slopes_aB = diff(aB_phi);
%         % check points where we're 5 times the mean, and large phi
%         new_remove = find((abs(slopes_aB'/mean(slopes_aB))>5).*(phi_vals(2:end)>2.6116));
% %         phi_vals(2:end)>2.6116
%         % also remove all points between the two limits
%         new_remove = min(new_remove):max(new_remove);
%         unphysical_removed = unique(cat(1,unphysical_removed, new_remove'));
        E_all_this(:,unphysical_removed) = nan;
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

    if clean_data
%         E = smooth(phi_vals, E, 'rloess');
%         [x,y] = prepareCurveData(phi_vals, E);
%         f = 
%         pp = csaps(phi_vals, E);
        [E,P] = csaps(phi_vals, E, 0.999999, phi_vals);
        if plot_curves
            figure(f1);
            plot(phi_vals, E);
%             P
%             0.999999
%             fnplt(pp);
        end
    end
    diff_E = diff(E);

%     if plot_curves
%         figure(f1);
%         hold on
%         plot(phi_vals,E,':')
%     end

    % find points where diff changes sign, these are stationary points.
    locs_pos = [];
    locs_neg = [];
    for ii = 2:length(diff_E)
        if ~(ignore_near_nan&&any(isnan(E(max([ii-2,1]):min([ii+2,length(E)])))))
            % maxima shouldn't really exist at small phi
            if (diff_E(ii)<0&&diff_E(ii-1)>0)&&ii>ceil(length(diff_E)/3)
                locs_pos = [locs_pos, ii];
            elseif diff_E(ii)>0&&diff_E(ii-1)<0
                % minima at small phi should be smaller than initial point
                if ii<ceil(length(diff_E)/2)
                    if E(ii)<E(1)
                        locs_neg = [locs_neg, ii];
                    end
                else
                    locs_neg = [locs_neg, ii];
                end
            end
        end
    end
    
    if plot_curves
        figure(f1);
        hold on
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
        E_near_min = E_all_this(1,min_range);
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
        [~,loc_neg_truly_min] = min(phi_vals(locs_neg));
        min_range =(locs_neg(loc_neg_truly_min)-1):(locs_neg(loc_neg_truly_min)+1);
        loc_min_index = locs_neg(loc_neg_truly_min);
        phi_near_min = phi_vals(min_range);
        E_near_min = E_all_this(1,min_range);
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





