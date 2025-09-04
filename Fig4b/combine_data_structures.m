function data_structure_combined = combine_data_structures(ds1, ds2)
    % firstly, check that the number of phi values is the same
    ab_size1 = size(ds1.alpha_B_vals_nonlinear);
    ab_size2 = size(ds2.alpha_B_vals_nonlinear);

    if ab_size1(1)~=ab_size2(1)
        error('number of phi vals not equal')
    end

    % get all fields from both
    fields1 = fieldnames(ds1);
    fields2 = fieldnames(ds2);

    % make sure they have the same fields
    if isempty(setdiff(fields1,fields2))
        % then loop over fields and combine
        for ii=1:length(fields1)
            field = fields1{ii};
            switch field
                case {'E_all_toroid','E_all_linear','E_all_nonlinear'}
                    data_structure_combined.(field) = cat(3,ds1.(field),ds2.(field));
                case {'parameter_set'}
                    data_structure_combined.(field) = cat(1,ds1.(field),ds2.(field));
                case {'alpha_A_vals_toroid','alpha_B_vals_toroid','S_A_vals_toroid',...
                        'S_B_vals_toroid', 'Sigma_vals_toroid', 'alpha_A_vals_linear',...
                        'alpha_B_vals_linear', 'S_A_vals_linear',  'S_B_vals_linear',...
                        'Sigma_vals_linear',  'alpha_A_vals_nonlinear', 'alpha_B_vals_nonlinear',...
                        'S_A_vals_nonlinear', 'S_B_vals_nonlinear', 'Sigma_vals_nonlinear',...
                        'rho_vals', 'phi_end_nonlinear'}
                    data_structure_combined.(field) = cat(2,ds1.(field),ds2.(field));
                case {'is_sorted', 'phi_vals'}
                    % don't do anything, they should be the same!
                    data_structure_combined.(field) = ds1.(field);
                otherwise
                    error('new field name, time to add it!')
            end
        end
    else
        error('fields are not the same between structures')
    end
end