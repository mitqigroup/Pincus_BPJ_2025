function name = get_curve_name(type)

curve_types = [1,2,3,4,5,6,0,7,8];
curve_names = ["C", "C_P", "P_C", "F_P", "F", "F_C", "0", "P", "8"];
d = dictionary(curve_types, curve_names);

name = d(type);