import os
import shutil

root = os.path.dirname(os.path.realpath(__file__))
line_to_change = "alpha_i_vals = alpha_i_vals_t(2);"
new_line = "alpha_i_vals = alpha_i_vals_t(1);"
# file_types =("nonlin_lookup.m", "sub_script.sh", "gridded_interp.mat")
file_types =("nonlin_lookup.m", "sub_script.sh")
dir_names = ("down", "middle")
oldline = "%phi_vals = [phi_vals(70:-1:1),phi_vals(71:length(phi_vals))];"
down_newline = "phi_vals = flip(phi_vals);"
middle_newline = "phi_vals = [phi_vals(70:-1:1),phi_vals(71:length(phi_vals))];"
base_name = "Sig_i_sweep"

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
    os.system('LLsub ./sub_script.sh [1,48,1] -T 24:00:00 -J "'+base_name+'_'+dir_name+'"')

# for root, dirs, files in os.walk(src_dir):
    # for file in files:
        # if file.endswith(file_types):
            # src_file = os.path.join(root, file)
            # dst_file = os.path.join(dst_dir, os.path.relpath(src_file, src_dir))
            # dst_folder = os.path.dirname(dst_file)
            # os.makedirs(dst_folder, exist_ok=True)
            # with open(src_file, 'r') as f:
                # lines = f.readlines()
            # with open(dst_file, 'w') as f:
                # for line in lines:
                    # if line.strip() == line_to_change:
                        # f.write(new_line + "\n")
                    # else:
                        # f.write(line)
            # if file.endswith(file_types[0]):
                # print(dst_folder[-1])
                # os.system('chmod u+x sub_script.sh')
                # os.system('LLsub ./sub_script.sh [1,48,1] -T 12:00:00 -J "case'+dst_folder[-1]+'_aa2rr2ss1"')
    # for d in dirs:
        # subdir = os.path.join(root, d)
        # os.chdir(subdir)
        # print(subdir)
        # print(d[-1])
