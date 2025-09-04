function outputs = get_all_global_mins(ds, xx, remove_points, slice, global_min_index)

yy = fit_spline(ds.phi_vals, ds.E_all_nonlinear(1,:,slice), xx, remove_points);
outputs.E_total_min = yy(global_min_index);

yy = fit_spline(ds.phi_vals, ds.E_all_nonlinear(2,:,slice), xx, remove_points);
outputs.E_adhesion_min = yy(global_min_index);

yy = fit_spline(ds.phi_vals, ds.E_all_nonlinear(3,:,slice), xx, remove_points);
outputs.E_stretchA_min = yy(global_min_index);

yy = fit_spline(ds.phi_vals, ds.E_all_nonlinear(4,:,slice), xx, remove_points);
outputs.E_stretchB_min = yy(global_min_index);

yy = fit_spline(ds.phi_vals, ds.E_all_nonlinear(5,:,slice), xx, remove_points);
outputs.E_bendA_min = yy(global_min_index);

yy = fit_spline(ds.phi_vals, ds.E_all_nonlinear(6,:,slice), xx, remove_points);
outputs.E_bendB_min = yy(global_min_index);

yy = fit_spline(ds.phi_vals, ds.alpha_A_vals_nonlinear(:,slice), xx, remove_points);
outputs.alpha_A_min = yy(global_min_index);

yy = fit_spline(ds.phi_vals, ds.alpha_B_vals_nonlinear(:,slice), xx, remove_points);
outputs.alpha_B_min = yy(global_min_index);

yy = fit_spline(ds.phi_vals, ds.S_A_vals_nonlinear(:,slice), xx, remove_points);
outputs.S_A_min = yy(global_min_index);

yy = fit_spline(ds.phi_vals, ds.S_B_vals_nonlinear(:,slice), xx, remove_points);
outputs.S_B_min = yy(global_min_index);

yy = fit_spline(ds.phi_vals, ds.Sigma_vals_nonlinear(:,slice), xx, remove_points);
outputs.Sigma_min = yy(global_min_index);














