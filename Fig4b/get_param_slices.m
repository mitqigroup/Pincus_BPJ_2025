function output_structure = get_param_slices(data_structure,slice)

size_out = size(data_structure.parameter_set);

output_structure.epsilon_slice = data_structure.parameter_set(slice,1);
output_structure.n0_slice = data_structure.parameter_set(slice,2);
output_structure.d_slice = data_structure.parameter_set(slice,3);
output_structure.R_slice = data_structure.parameter_set(slice,4);
output_structure.kD_slice = data_structure.parameter_set(slice,5);
output_structure.kappa_slice = data_structure.parameter_set(slice,6);
output_structure.alpha_i_slice = data_structure.parameter_set(slice,7);
if size_out(2)==8
    output_structure.m_slice = data_structure.parameter_set(slice,8);
end
output_structure.sigma_slice = pi*output_structure.R_slice.^2./output_structure.d_slice.^2;
output_structure.zeta_slice = -output_structure.epsilon_slice.*output_structure.n0_slice./output_structure.kD_slice;
output_structure.Sig_bar_init_slice = output_structure.kD_slice.*output_structure.alpha_i_slice.*...
    output_structure.R_slice.^2./output_structure.kappa_slice;
output_structure.w_bar_slice = -2*output_structure.epsilon_slice.*output_structure.n0_slice.*...
    output_structure.R_slice.^2./output_structure.kappa_slice;
output_structure.R_lam_k_slice = output_structure.R_slice.*...
    sqrt(output_structure.kappa_slice./output_structure.kD_slice);

end