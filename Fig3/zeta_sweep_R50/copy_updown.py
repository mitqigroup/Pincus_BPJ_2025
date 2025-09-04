import os
import shutil

root = "/home/gridsan/ipincus/membranes/nonlinear_parameter_sweeps/ep_sig_sweep/zeta_sweep_R50"
# line_to_change = "alpha_i_vals = alpha_i_vals_t(2);"
# new_line = "alpha_i_vals = alpha_i_vals_t(1);"
# file_types =("nonlin_lookup.m", "sub_script.sh", "gridded_interp.mat")
file_types =("nonlin_lookup.m", "sub_script.sh")
dir_names = ("down", "middle")
oldline = "%phi_vals = [phi_vals(70:-1:1),phi_vals(71:length(phi_vals))];"
down_newline = "phi_vals = flip(phi_vals);"
middle_newline = "phi_vals = [phi_vals(70:-1:1),phi_vals(71:length(phi_vals))];"
base_name = "Comprehensive_test_updown"

for dir_name in dir_names:
    dir_full = os.path.join(root, dir_name)
    os.makedirs(dir_full, exist_ok=True)
    os.chdir(dir_full)
    for file in file_types:
        file_src_path = os.path.join(root, file)
        file_dst_path = os.path.join(dir_full, file)
        shutil.copy2(file_src_path, file_dst_path)
        if os.path.basename(file_dst_path)=="nonlin_lookup.m" and dir_name=="down":
            with open(file_dst_path, "r") as dst_file:
                lines = dst_file.readlines()
            with open(file_dst_path, "w") as dst_file:
                for line in lines:
                    if line.strip() == oldline:
                        dst_file.write(down_newline + "\n")
                    else:
                        dst_file.write(line)
        elif os.path.basename(file_dst_path)=="nonlin_lookup.m" and dir_name=="middle":
            with open(file_dst_path, "r") as dst_file:
                lines = dst_file.readlines()
            with open(file_dst_path, "w") as dst_file:
                for line in lines:
                    if line.strip() == oldline:
                        dst_file.write(middle_newline + "\n")
                    else:
                        dst_file.write(line)
    os.system('LLsub ./sub_script.sh [1,48,1] -T 4:00:00 -J "'+base_name+'_'+dir_name+'"')
