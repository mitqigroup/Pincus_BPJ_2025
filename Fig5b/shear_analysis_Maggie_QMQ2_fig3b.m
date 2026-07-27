clear variables

[filePath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
current_location = [filePath,'/middle_QMQ17'];
ds_middle_1 = read_and_unpermute(current_location, 1, 'sortvalues', true);
%ds_middle = read_and_unpermute(current_location, 96, 'sortvalues', true);
current_location = [filePath,'/down_QMQ17'];
ds_down_1 = read_and_unpermute(current_location, 1, 'sortvalues', true);
%ds_down = read_and_unpermute(current_location, 96, 'sortvalues', true);

current_location = [filePath,'/middle_QMQ17_long'];
ds_middle_2 = read_and_unpermute(current_location, 1, 'sortvalues', true);
current_location = [filePath,'/down_QMQ17_long'];
ds_down_2 = read_and_unpermute(current_location, 1, 'sortvalues', true);

current_location = [filePath,'/middle_QMQ17_long2'];
ds_middle_3 = read_and_unpermute(current_location, 1, 'sortvalues', true);
current_location = [filePath,'/down_QMQ17_long2'];
ds_down_3 = read_and_unpermute(current_location, 1, 'sortvalues', true);
newcolors = flip(["#003f5c", "#374c80", "#7a5195", "#bc5090", "#ef5675", ...
    "#ff764a", "#ffa600"]);

current_location = [filePath,'/middle_QMQ17_long3'];
ds_middle_4 = read_and_unpermute(current_location, 1, 'sortvalues', true);
current_location = [filePath,'/down_QMQ17_long3'];
ds_down_4 = read_and_unpermute(current_location, 1, 'sortvalues', true);

current_location = [filePath,'/middle_QMQ17_long4'];
ds_middle_5 = read_and_unpermute(current_location, 1, 'sortvalues', true);
current_location = [filePath,'/down_QMQ17_long4'];
ds_down_5 = read_and_unpermute(current_location, 1, 'sortvalues', true);

current_location = [filePath,'/middle_QMQ17_long5'];
ds_middle_6 = read_and_unpermute(current_location, 1, 'sortvalues', true);
current_location = [filePath,'/down_QMQ17_long5'];
ds_down_6 = read_and_unpermute(current_location, 1, 'sortvalues', true);

ds_middle=ds_middle_1;
ds_down=ds_down_1;
ds_middle.parameter_set=zeros(size( ds_middle_1.parameter_set,1)*4+18,size( ds_middle_1.parameter_set,2));
ds_down.parameter_set=zeros(size( ds_down_1.parameter_set,1)*4+18,size( ds_down_1.parameter_set,2));
ds_middle.E_all_toroid=zeros(size( ds_middle_1.E_all_toroid,1),size( ds_middle_1.E_all_toroid,2),4*size( ds_middle_1.E_all_toroid,3)+18);
ds_down.E_all_toroid=zeros(size( ds_middle_1.E_all_toroid,1),size( ds_middle_1.E_all_toroid,2),4*size( ds_middle_1.E_all_toroid,3)+18);
ds_middle.alpha_A_vals_toroid=zeros(size( ds_middle_1.alpha_A_vals_toroid,1),size( ds_middle_1.alpha_A_vals_toroid,2)*4+18);
ds_down.alpha_A_vals_toroid=zeros(size( ds_down_1.alpha_A_vals_toroid,1),size( ds_down_1.alpha_A_vals_toroid,2)*4+18);
ds_middle.alpha_B_vals_toroid=zeros(size( ds_middle_1.alpha_B_vals_toroid,1),size( ds_middle_1.alpha_B_vals_toroid,2)*4+18);
ds_down.alpha_B_vals_toroid=zeros(size( ds_down_1.alpha_B_vals_toroid,1),size( ds_down_1.alpha_B_vals_toroid,2)*4+18);
ds_middle.S_A_vals_toroid=zeros(size( ds_middle_1.S_A_vals_toroid,1),size( ds_middle_1.S_A_vals_toroid,2)*4+18);
ds_down.S_A_vals_toroid=zeros(size( ds_down_1.S_A_vals_toroid,1),size( ds_down_1.S_A_vals_toroid,2)*4+18);
ds_middle.S_B_vals_toroid=zeros(size( ds_middle_1.S_B_vals_toroid,1),size( ds_middle_1.S_B_vals_toroid,2)*4+18);
ds_down.S_B_vals_toroid=zeros(size( ds_down_1.S_B_vals_toroid,1),size( ds_down_1.S_B_vals_toroid,2)*4+18);
ds_middle.Sigma_vals_toroid=zeros(size( ds_middle_1.Sigma_vals_toroid,1),size( ds_middle_1.Sigma_vals_toroid,2)*4+18);
ds_down.Sigma_vals_toroid=zeros(size( ds_down_1.Sigma_vals_toroid,1),size( ds_down_1.Sigma_vals_toroid,2)*4+18);
ds_middle.rho_vals=zeros(size( ds_middle_1.rho_vals,1),size( ds_middle_1.rho_vals,2)*4+18);
ds_down.rho_vals=zeros(size( ds_down_1.rho_vals,1),size( ds_down_1.rho_vals,2)*4+18);
ds_middle.E_all_linear=zeros(size( ds_middle_1.E_all_linear,1),size( ds_middle_1.E_all_linear,2),4*size( ds_middle_1.E_all_linear,3)+18);
ds_down.E_all_linear=zeros(size( ds_middle_1.E_all_linear,1),size( ds_middle_1.E_all_linear,2),4*size( ds_middle_1.E_all_linear,3)+18);
ds_middle.alpha_A_vals_linear=zeros(size( ds_middle_1.alpha_A_vals_linear,1),size( ds_middle_1.alpha_A_vals_linear,2)*4+18);
ds_down.alpha_A_vals_linear=zeros(size( ds_down_1.alpha_A_vals_linear,1),size( ds_down_1.alpha_A_vals_linear,2)*4+18);
ds_middle.alpha_B_vals_linear=zeros(size( ds_middle_1.alpha_B_vals_linear,1),size( ds_middle_1.alpha_B_vals_linear,2)*4+18);
ds_down.alpha_B_vals_linear=zeros(size( ds_down_1.alpha_B_vals_linear,1),size( ds_down_1.alpha_B_vals_linear,2)*4+18);
ds_middle.S_A_vals_linear=zeros(size( ds_middle_1.S_A_vals_linear,1),size( ds_middle_1.S_A_vals_linear,2)*4+18);
ds_down.S_A_vals_linear=zeros(size( ds_down_1.S_A_vals_linear,1),size( ds_down_1.S_A_vals_linear,2)*4+18);
ds_middle.S_B_vals_linear=zeros(size( ds_middle_1.S_B_vals_linear,1),size( ds_middle_1.S_B_vals_linear,2)*4+18);
ds_down.S_B_vals_linear=zeros(size( ds_down_1.S_B_vals_linear,1),size( ds_down_1.S_B_vals_linear,2)*4+18);
ds_middle.Sigma_vals_linear=zeros(size( ds_middle_1.Sigma_vals_linear,1),size( ds_middle_1.Sigma_vals_linear,2)*4+18);
ds_down.Sigma_vals_linear=zeros(size( ds_down_1.Sigma_vals_linear,1),size( ds_down_1.Sigma_vals_linear,2)*4+18);
ds_middle.E_all_nonlinear=zeros(size( ds_middle_1.E_all_nonlinear,1),size( ds_middle_1.E_all_nonlinear,2),4*size( ds_middle_1.E_all_nonlinear,3)+18);
ds_down.E_all_nonlinear=zeros(size( ds_middle_1.E_all_nonlinear,1),size( ds_middle_1.E_all_nonlinear,2),4*size( ds_middle_1.E_all_nonlinear,3)+18);
ds_middle.alpha_A_vals_nonlinear=zeros(size( ds_middle_1.alpha_A_vals_nonlinear,1),size( ds_middle_1.alpha_A_vals_nonlinear,2)*4+18);
ds_down.alpha_A_vals_nonlinear=zeros(size( ds_down_1.alpha_A_vals_nonlinear,1),size( ds_down_1.alpha_A_vals_nonlinear,2)*4+18);
ds_middle.alpha_B_vals_nonlinear=zeros(size( ds_middle_1.alpha_B_vals_nonlinear,1),size( ds_middle_1.alpha_B_vals_nonlinear,2)*4+18);
ds_down.alpha_B_vals_nonlinear=zeros(size( ds_down_1.alpha_B_vals_nonlinear,1),size( ds_down_1.alpha_B_vals_nonlinear,2)*4+18);
ds_middle.S_A_vals_nonlinear=zeros(size( ds_middle_1.S_A_vals_nonlinear,1),size( ds_middle_1.S_A_vals_nonlinear,2)*4+18);
ds_down.S_A_vals_nonlinear=zeros(size( ds_down_1.S_A_vals_nonlinear,1),size( ds_down_1.S_A_vals_nonlinear,2)*4+18);
ds_middle.S_B_vals_nonlinear=zeros(size( ds_middle_1.S_B_vals_nonlinear,1),size( ds_middle_1.S_B_vals_nonlinear,2)*4+18);
ds_down.S_B_vals_nonlinear=zeros(size( ds_down_1.S_B_vals_nonlinear,1),size( ds_down_1.S_B_vals_nonlinear,2)*4+18);
ds_middle.Sigma_vals_nonlinear=zeros(size( ds_middle_1.Sigma_vals_nonlinear,1),size( ds_middle_1.Sigma_vals_nonlinear,2)*4+18);
ds_down.Sigma_vals_nonlinear=zeros(size( ds_down_1.Sigma_vals_nonlinear,1),size( ds_down_1.Sigma_vals_nonlinear,2)*4+18);
ds_middle.phi_end_nonlinear=zeros(size( ds_middle_1.phi_end_nonlinear,1),size( ds_middle_1.phi_end_nonlinear,2)*4+18);
ds_down.phi_end_nonlinear=zeros(size( ds_down_1.phi_end_nonlinear,1),size( ds_down_1.phi_end_nonlinear,2)*4+18);


R_length=7;
sigma_length=4;
kD_length=1;
epsilon_length=1;
n0_length=1;
kappa_length=1;
alpha_i_length=1;

ii=0;
ii_1=0;
iivec_1=zeros(size( ds_middle.parameter_set,1),1);
ii_2=0;
iivec_2=zeros(size( ds_middle.parameter_set,1),1);
ii_3=0;
ii_4=0;
ii_5=0;
ii_6=0;

sortedRvec=sort([1:5:68 13:5:53]);
for rr = 1:length(sortedRvec)
    for ss = 1:2*sigma_length
        for kk = 1:kD_length
            for ee = 1:epsilon_length
                for nn = 1:n0_length
                    for pp = 1:kappa_length
                        for aa = 1:alpha_i_length
                            ii = ii+1;
                            if rem(sortedRvec(rr)-3,5)==0&&rem(ss,2)==1
                                ii_5=ii_5+1;
                                ds_middle.parameter_set(ii,:)=ds_middle_5.parameter_set(ii_5,:);
                                ds_down.parameter_set(ii,:)=ds_down_5.parameter_set(ii_5,:);
                                ds_middle.E_all_toroid(:,:,ii)=ds_middle_5.E_all_toroid(:,:,ii_5);
                                ds_down.E_all_toroid(:,:,ii)=ds_down_5.E_all_toroid(:,:,ii_5);
                                ds_middle.alpha_A_vals_toroid(:,ii)=ds_middle_5.alpha_A_vals_toroid(:,ii_5);
                                ds_down.alpha_A_vals_toroid(:,ii)=ds_down_5.alpha_A_vals_toroid(:,ii_5);
                                ds_middle.alpha_B_vals_toroid(:,ii)=ds_middle_5.alpha_B_vals_toroid(:,ii_5);
                                ds_down.alpha_B_vals_toroid(:,ii)=ds_down_5.alpha_B_vals_toroid(:,ii_5);
                                ds_middle.S_A_vals_toroid(:,ii)=ds_middle_5.S_A_vals_toroid(:,ii_5);
                                ds_down.S_A_vals_toroid(:,ii)=ds_down_5.S_A_vals_toroid(:,ii_5);
                                ds_middle.S_B_vals_toroid(:,ii)=ds_middle_5.S_B_vals_toroid(:,ii_5);
                                ds_down.S_B_vals_toroid(:,ii)=ds_down_5.S_B_vals_toroid(:,ii_5);
                                ds_middle.Sigma_vals_toroid(:,ii)=ds_middle_5.Sigma_vals_toroid(:,ii_5);
                                ds_down.Sigma_vals_toroid(:,ii)=ds_down_5.Sigma_vals_toroid(:,ii_5);
                                ds_middle.rho_vals(:,ii)=ds_middle_5.rho_vals(:,ii_5);
                                ds_down.rho_vals(:,ii)=ds_down_5.rho_vals(:,ii_5);
                                ds_middle.E_all_linear(:,:,ii)=ds_middle_5.E_all_linear(:,:,ii_5);
                                ds_down.E_all_linear(:,:,ii)=ds_down_5.E_all_linear(:,:,ii_5);
                                ds_middle.alpha_A_vals_linear(:,ii)=ds_middle_5.alpha_A_vals_linear(:,ii_5);
                                ds_down.alpha_A_vals_linear(:,ii)=ds_down_5.alpha_A_vals_linear(:,ii_5);
                                ds_middle.alpha_B_vals_linear(:,ii)=ds_middle_5.alpha_B_vals_linear(:,ii_5);
                                ds_down.alpha_B_vals_linear(:,ii)=ds_down_5.alpha_B_vals_linear(:,ii_5);
                                ds_middle.S_A_vals_linear(:,ii)=ds_middle_5.S_A_vals_linear(:,ii_5);
                                ds_down.S_A_vals_linear(:,ii)=ds_down_5.S_A_vals_linear(:,ii_5);
                                ds_middle.S_B_vals_linear(:,ii)=ds_middle_5.S_B_vals_linear(:,ii_5);
                                ds_down.S_B_vals_linear(:,ii)=ds_down_5.S_B_vals_linear(:,ii_5);
                                ds_middle.Sigma_vals_linear(:,ii)=ds_middle_5.Sigma_vals_linear(:,ii_5);
                                ds_down.Sigma_vals_linear(:,ii)=ds_down_5.Sigma_vals_linear(:,ii_5);
                                ds_middle.E_all_nonlinear(:,:,ii)=ds_middle_5.E_all_nonlinear(:,:,ii_5);
                                ds_down.E_all_nonlinear(:,:,ii)=ds_down_5.E_all_nonlinear(:,:,ii_5);
                                ds_middle.alpha_A_vals_nonlinear(:,ii)=ds_middle_5.alpha_A_vals_nonlinear(:,ii_5);
                                ds_down.alpha_A_vals_nonlinear(:,ii)=ds_down_5.alpha_A_vals_nonlinear(:,ii_5);
                                ds_middle.alpha_B_vals_nonlinear(:,ii)=ds_middle_5.alpha_B_vals_nonlinear(:,ii_5);
                                ds_down.alpha_B_vals_nonlinear(:,ii)=ds_down_5.alpha_B_vals_nonlinear(:,ii_5);
                                ds_middle.S_A_vals_nonlinear(:,ii)=ds_middle_5.S_A_vals_nonlinear(:,ii_5);
                                ds_down.S_A_vals_nonlinear(:,ii)=ds_down_5.S_A_vals_nonlinear(:,ii_5);
                                ds_middle.S_B_vals_nonlinear(:,ii)=ds_middle_5.S_B_vals_nonlinear(:,ii_5);
                                ds_down.S_B_vals_nonlinear(:,ii)=ds_down_5.S_B_vals_nonlinear(:,ii_5);
                                ds_middle.Sigma_vals_nonlinear(:,ii)=ds_middle_5.Sigma_vals_nonlinear(:,ii_5);
                                ds_down.Sigma_vals_nonlinear(:,ii)=ds_down_5.Sigma_vals_nonlinear(:,ii_5);
                                ds_middle.phi_end_nonlinear(:,ii)=ds_middle_5.phi_end_nonlinear(:,ii_5);
                                ds_down.phi_end_nonlinear(:,ii)=ds_down_5.phi_end_nonlinear(:,ii_5);
                            elseif rem(sortedRvec(rr)-3,5)==0&&rem(ss,2)==0
                                ii_6=ii_6+1;
                                ds_middle.parameter_set(ii,:)=ds_middle_6.parameter_set(ii_6,:);
                                ds_down.parameter_set(ii,:)=ds_down_6.parameter_set(ii_6,:);
                                ds_middle.E_all_toroid(:,:,ii)=ds_middle_6.E_all_toroid(:,:,ii_6);
                                ds_down.E_all_toroid(:,:,ii)=ds_down_6.E_all_toroid(:,:,ii_6);
                                ds_middle.alpha_A_vals_toroid(:,ii)=ds_middle_6.alpha_A_vals_toroid(:,ii_6);
                                ds_down.alpha_A_vals_toroid(:,ii)=ds_down_6.alpha_A_vals_toroid(:,ii_6);
                                ds_middle.alpha_B_vals_toroid(:,ii)=ds_middle_6.alpha_B_vals_toroid(:,ii_6);
                                ds_down.alpha_B_vals_toroid(:,ii)=ds_down_6.alpha_B_vals_toroid(:,ii_6);
                                ds_middle.S_A_vals_toroid(:,ii)=ds_middle_6.S_A_vals_toroid(:,ii_6);
                                ds_down.S_A_vals_toroid(:,ii)=ds_down_6.S_A_vals_toroid(:,ii_6);
                                ds_middle.S_B_vals_toroid(:,ii)=ds_middle_6.S_B_vals_toroid(:,ii_6);
                                ds_down.S_B_vals_toroid(:,ii)=ds_down_6.S_B_vals_toroid(:,ii_6);
                                ds_middle.Sigma_vals_toroid(:,ii)=ds_middle_6.Sigma_vals_toroid(:,ii_6);
                                ds_down.Sigma_vals_toroid(:,ii)=ds_down_6.Sigma_vals_toroid(:,ii_6);
                                ds_middle.rho_vals(:,ii)=ds_middle_6.rho_vals(:,ii_6);
                                ds_down.rho_vals(:,ii)=ds_down_6.rho_vals(:,ii_6);
                                ds_middle.E_all_linear(:,:,ii)=ds_middle_6.E_all_linear(:,:,ii_6);
                                ds_down.E_all_linear(:,:,ii)=ds_down_6.E_all_linear(:,:,ii_6);
                                ds_middle.alpha_A_vals_linear(:,ii)=ds_middle_6.alpha_A_vals_linear(:,ii_6);
                                ds_down.alpha_A_vals_linear(:,ii)=ds_down_6.alpha_A_vals_linear(:,ii_6);
                                ds_middle.alpha_B_vals_linear(:,ii)=ds_middle_6.alpha_B_vals_linear(:,ii_6);
                                ds_down.alpha_B_vals_linear(:,ii)=ds_down_6.alpha_B_vals_linear(:,ii_6);
                                ds_middle.S_A_vals_linear(:,ii)=ds_middle_6.S_A_vals_linear(:,ii_6);
                                ds_down.S_A_vals_linear(:,ii)=ds_down_6.S_A_vals_linear(:,ii_6);
                                ds_middle.S_B_vals_linear(:,ii)=ds_middle_6.S_B_vals_linear(:,ii_6);
                                ds_down.S_B_vals_linear(:,ii)=ds_down_6.S_B_vals_linear(:,ii_6);
                                ds_middle.Sigma_vals_linear(:,ii)=ds_middle_6.Sigma_vals_linear(:,ii_6);
                                ds_down.Sigma_vals_linear(:,ii)=ds_down_6.Sigma_vals_linear(:,ii_6);
                                ds_middle.E_all_nonlinear(:,:,ii)=ds_middle_6.E_all_nonlinear(:,:,ii_6);
                                ds_down.E_all_nonlinear(:,:,ii)=ds_down_6.E_all_nonlinear(:,:,ii_6);
                                ds_middle.alpha_A_vals_nonlinear(:,ii)=ds_middle_6.alpha_A_vals_nonlinear(:,ii_6);
                                ds_down.alpha_A_vals_nonlinear(:,ii)=ds_down_6.alpha_A_vals_nonlinear(:,ii_6);
                                ds_middle.alpha_B_vals_nonlinear(:,ii)=ds_middle_6.alpha_B_vals_nonlinear(:,ii_6);
                                ds_down.alpha_B_vals_nonlinear(:,ii)=ds_down_6.alpha_B_vals_nonlinear(:,ii_6);
                                ds_middle.S_A_vals_nonlinear(:,ii)=ds_middle_6.S_A_vals_nonlinear(:,ii_6);
                                ds_down.S_A_vals_nonlinear(:,ii)=ds_down_6.S_A_vals_nonlinear(:,ii_6);
                                ds_middle.S_B_vals_nonlinear(:,ii)=ds_middle_6.S_B_vals_nonlinear(:,ii_6);
                                ds_down.S_B_vals_nonlinear(:,ii)=ds_down_6.S_B_vals_nonlinear(:,ii_6);
                                ds_middle.Sigma_vals_nonlinear(:,ii)=ds_middle_6.Sigma_vals_nonlinear(:,ii_6);
                                ds_down.Sigma_vals_nonlinear(:,ii)=ds_down_6.Sigma_vals_nonlinear(:,ii_6);
                                ds_middle.phi_end_nonlinear(:,ii)=ds_middle_6.phi_end_nonlinear(:,ii_6);
                                ds_down.phi_end_nonlinear(:,ii)=ds_down_6.phi_end_nonlinear(:,ii_6);
                            elseif rem(sortedRvec(rr),2)==1&&rem(ss,2)==1
                                ii_1=ii_1+1;
                                ds_middle.parameter_set(ii,:)=ds_middle_1.parameter_set(ii_1,:);
                                ds_down.parameter_set(ii,:)=ds_down_1.parameter_set(ii_1,:);
                                ds_middle.E_all_toroid(:,:,ii)=ds_middle_1.E_all_toroid(:,:,ii_1);
                                ds_down.E_all_toroid(:,:,ii)=ds_down_1.E_all_toroid(:,:,ii_1);
                                ds_middle.alpha_A_vals_toroid(:,ii)=ds_middle_1.alpha_A_vals_toroid(:,ii_1);
                                ds_down.alpha_A_vals_toroid(:,ii)=ds_down_1.alpha_A_vals_toroid(:,ii_1);
                                ds_middle.alpha_B_vals_toroid(:,ii)=ds_middle_1.alpha_B_vals_toroid(:,ii_1);
                                ds_down.alpha_B_vals_toroid(:,ii)=ds_down_1.alpha_B_vals_toroid(:,ii_1);
                                ds_middle.S_A_vals_toroid(:,ii)=ds_middle_1.S_A_vals_toroid(:,ii_1);
                                ds_down.S_A_vals_toroid(:,ii)=ds_down_1.S_A_vals_toroid(:,ii_1);
                                ds_middle.S_B_vals_toroid(:,ii)=ds_middle_1.S_B_vals_toroid(:,ii_1);
                                ds_down.S_B_vals_toroid(:,ii)=ds_down_1.S_B_vals_toroid(:,ii_1);
                                ds_middle.Sigma_vals_toroid(:,ii)=ds_middle_1.Sigma_vals_toroid(:,ii_1);
                                ds_down.Sigma_vals_toroid(:,ii)=ds_down_1.Sigma_vals_toroid(:,ii_1);
                                ds_middle.rho_vals(:,ii)=ds_middle_1.rho_vals(:,ii_1);
                                ds_down.rho_vals(:,ii)=ds_down_1.rho_vals(:,ii_1);
                                ds_middle.E_all_linear(:,:,ii)=ds_middle_1.E_all_linear(:,:,ii_1);
                                ds_down.E_all_linear(:,:,ii)=ds_down_1.E_all_linear(:,:,ii_1);
                                ds_middle.alpha_A_vals_linear(:,ii)=ds_middle_1.alpha_A_vals_linear(:,ii_1);
                                ds_down.alpha_A_vals_linear(:,ii)=ds_down_1.alpha_A_vals_linear(:,ii_1);
                                ds_middle.alpha_B_vals_linear(:,ii)=ds_middle_1.alpha_B_vals_linear(:,ii_1);
                                ds_down.alpha_B_vals_linear(:,ii)=ds_down_1.alpha_B_vals_linear(:,ii_1);
                                ds_middle.S_A_vals_linear(:,ii)=ds_middle_1.S_A_vals_linear(:,ii_1);
                                ds_down.S_A_vals_linear(:,ii)=ds_down_1.S_A_vals_linear(:,ii_1);
                                ds_middle.S_B_vals_linear(:,ii)=ds_middle_1.S_B_vals_linear(:,ii_1);
                                ds_down.S_B_vals_linear(:,ii)=ds_down_1.S_B_vals_linear(:,ii_1);
                                ds_middle.Sigma_vals_linear(:,ii)=ds_middle_1.Sigma_vals_linear(:,ii_1);
                                ds_down.Sigma_vals_linear(:,ii)=ds_down_1.Sigma_vals_linear(:,ii_1);
                                ds_middle.E_all_nonlinear(:,:,ii)=ds_middle_1.E_all_nonlinear(:,:,ii_1);
                                ds_down.E_all_nonlinear(:,:,ii)=ds_down_1.E_all_nonlinear(:,:,ii_1);
                                ds_middle.alpha_A_vals_nonlinear(:,ii)=ds_middle_1.alpha_A_vals_nonlinear(:,ii_1);
                                ds_down.alpha_A_vals_nonlinear(:,ii)=ds_down_1.alpha_A_vals_nonlinear(:,ii_1);
                                ds_middle.alpha_B_vals_nonlinear(:,ii)=ds_middle_1.alpha_B_vals_nonlinear(:,ii_1);
                                ds_down.alpha_B_vals_nonlinear(:,ii)=ds_down_1.alpha_B_vals_nonlinear(:,ii_1);
                                ds_middle.S_A_vals_nonlinear(:,ii)=ds_middle_1.S_A_vals_nonlinear(:,ii_1);
                                ds_down.S_A_vals_nonlinear(:,ii)=ds_down_1.S_A_vals_nonlinear(:,ii_1);
                                ds_middle.S_B_vals_nonlinear(:,ii)=ds_middle_1.S_B_vals_nonlinear(:,ii_1);
                                ds_down.S_B_vals_nonlinear(:,ii)=ds_down_1.S_B_vals_nonlinear(:,ii_1);
                                ds_middle.Sigma_vals_nonlinear(:,ii)=ds_middle_1.Sigma_vals_nonlinear(:,ii_1);
                                ds_down.Sigma_vals_nonlinear(:,ii)=ds_down_1.Sigma_vals_nonlinear(:,ii_1);
                                ds_middle.phi_end_nonlinear(:,ii)=ds_middle_1.phi_end_nonlinear(:,ii_1);
                                ds_down.phi_end_nonlinear(:,ii)=ds_down_1.phi_end_nonlinear(:,ii_1);
                            elseif rem(sortedRvec(rr),2)==1&&rem(ss,2)==0
                                ii_2=ii_2+1;
                                ds_middle.parameter_set(ii,:)=ds_middle_2.parameter_set(ii_2,:);
                                ds_down.parameter_set(ii,:)=ds_down_2.parameter_set(ii_2,:);
                                ds_middle.E_all_toroid(:,:,ii)=ds_middle_2.E_all_toroid(:,:,ii_2);
                                ds_down.E_all_toroid(:,:,ii)=ds_down_2.E_all_toroid(:,:,ii_2);
                                ds_middle.alpha_A_vals_toroid(:,ii)=ds_middle_2.alpha_A_vals_toroid(:,ii_2);
                                ds_down.alpha_A_vals_toroid(:,ii)=ds_down_2.alpha_A_vals_toroid(:,ii_2);
                                ds_middle.alpha_B_vals_toroid(:,ii)=ds_middle_2.alpha_B_vals_toroid(:,ii_2);
                                ds_down.alpha_B_vals_toroid(:,ii)=ds_down_2.alpha_B_vals_toroid(:,ii_2);
                                ds_middle.S_A_vals_toroid(:,ii)=ds_middle_2.S_A_vals_toroid(:,ii_2);
                                ds_down.S_A_vals_toroid(:,ii)=ds_down_2.S_A_vals_toroid(:,ii_2);
                                ds_middle.S_B_vals_toroid(:,ii)=ds_middle_2.S_B_vals_toroid(:,ii_2);
                                ds_down.S_B_vals_toroid(:,ii)=ds_down_2.S_B_vals_toroid(:,ii_2);
                                ds_middle.Sigma_vals_toroid(:,ii)=ds_middle_2.Sigma_vals_toroid(:,ii_2);
                                ds_down.Sigma_vals_toroid(:,ii)=ds_down_2.Sigma_vals_toroid(:,ii_2);
                                ds_middle.rho_vals(:,ii)=ds_middle_2.rho_vals(:,ii_2);
                                ds_down.rho_vals(:,ii)=ds_down_2.rho_vals(:,ii_2);
                                ds_middle.E_all_linear(:,:,ii)=ds_middle_2.E_all_linear(:,:,ii_2);
                                ds_down.E_all_linear(:,:,ii)=ds_down_2.E_all_linear(:,:,ii_2);
                                ds_middle.alpha_A_vals_linear(:,ii)=ds_middle_2.alpha_A_vals_linear(:,ii_2);
                                ds_down.alpha_A_vals_linear(:,ii)=ds_down_2.alpha_A_vals_linear(:,ii_2);
                                ds_middle.alpha_B_vals_linear(:,ii)=ds_middle_2.alpha_B_vals_linear(:,ii_2);
                                ds_down.alpha_B_vals_linear(:,ii)=ds_down_2.alpha_B_vals_linear(:,ii_2);
                                ds_middle.S_A_vals_linear(:,ii)=ds_middle_2.S_A_vals_linear(:,ii_2);
                                ds_down.S_A_vals_linear(:,ii)=ds_down_2.S_A_vals_linear(:,ii_2);
                                ds_middle.S_B_vals_linear(:,ii)=ds_middle_2.S_B_vals_linear(:,ii_2);
                                ds_down.S_B_vals_linear(:,ii)=ds_down_2.S_B_vals_linear(:,ii_2);
                                ds_middle.Sigma_vals_linear(:,ii)=ds_middle_2.Sigma_vals_linear(:,ii_2);
                                ds_down.Sigma_vals_linear(:,ii)=ds_down_2.Sigma_vals_linear(:,ii_2);
                                ds_middle.E_all_nonlinear(:,:,ii)=ds_middle_2.E_all_nonlinear(:,:,ii_2);
                                ds_down.E_all_nonlinear(:,:,ii)=ds_down_2.E_all_nonlinear(:,:,ii_2);
                                ds_middle.alpha_A_vals_nonlinear(:,ii)=ds_middle_2.alpha_A_vals_nonlinear(:,ii_2);
                                ds_down.alpha_A_vals_nonlinear(:,ii)=ds_down_2.alpha_A_vals_nonlinear(:,ii_2);
                                ds_middle.alpha_B_vals_nonlinear(:,ii)=ds_middle_2.alpha_B_vals_nonlinear(:,ii_2);
                                ds_down.alpha_B_vals_nonlinear(:,ii)=ds_down_2.alpha_B_vals_nonlinear(:,ii_2);
                                ds_middle.S_A_vals_nonlinear(:,ii)=ds_middle_2.S_A_vals_nonlinear(:,ii_2);
                                ds_down.S_A_vals_nonlinear(:,ii)=ds_down_2.S_A_vals_nonlinear(:,ii_2);
                                ds_middle.S_B_vals_nonlinear(:,ii)=ds_middle_2.S_B_vals_nonlinear(:,ii_2);
                                ds_down.S_B_vals_nonlinear(:,ii)=ds_down_2.S_B_vals_nonlinear(:,ii_2);
                                ds_middle.Sigma_vals_nonlinear(:,ii)=ds_middle_2.Sigma_vals_nonlinear(:,ii_2);
                                ds_down.Sigma_vals_nonlinear(:,ii)=ds_down_2.Sigma_vals_nonlinear(:,ii_2);
                                ds_middle.phi_end_nonlinear(:,ii)=ds_middle_2.phi_end_nonlinear(:,ii_2);
                                ds_down.phi_end_nonlinear(:,ii)=ds_down_2.phi_end_nonlinear(:,ii_2);
                            elseif rem(sortedRvec(rr),2)==0&&rem(ss,2)==1
                                ii_3=ii_3+1;
                                ds_middle.parameter_set(ii,:)=ds_middle_3.parameter_set(ii_3,:);
                                ds_down.parameter_set(ii,:)=ds_down_3.parameter_set(ii_3,:);
                                ds_middle.E_all_toroid(:,:,ii)=ds_middle_3.E_all_toroid(:,:,ii_3);
                                ds_down.E_all_toroid(:,:,ii)=ds_down_3.E_all_toroid(:,:,ii_3);
                                ds_middle.alpha_A_vals_toroid(:,ii)=ds_middle_3.alpha_A_vals_toroid(:,ii_3);
                                ds_down.alpha_A_vals_toroid(:,ii)=ds_down_3.alpha_A_vals_toroid(:,ii_3);
                                ds_middle.alpha_B_vals_toroid(:,ii)=ds_middle_3.alpha_B_vals_toroid(:,ii_3);
                                ds_down.alpha_B_vals_toroid(:,ii)=ds_down_3.alpha_B_vals_toroid(:,ii_3);
                                ds_middle.S_A_vals_toroid(:,ii)=ds_middle_3.S_A_vals_toroid(:,ii_3);
                                ds_down.S_A_vals_toroid(:,ii)=ds_down_3.S_A_vals_toroid(:,ii_3);
                                ds_middle.S_B_vals_toroid(:,ii)=ds_middle_3.S_B_vals_toroid(:,ii_3);
                                ds_down.S_B_vals_toroid(:,ii)=ds_down_3.S_B_vals_toroid(:,ii_3);
                                ds_middle.Sigma_vals_toroid(:,ii)=ds_middle_3.Sigma_vals_toroid(:,ii_3);
                                ds_down.Sigma_vals_toroid(:,ii)=ds_down_3.Sigma_vals_toroid(:,ii_3);
                                ds_middle.rho_vals(:,ii)=ds_middle_3.rho_vals(:,ii_3);
                                ds_down.rho_vals(:,ii)=ds_down_3.rho_vals(:,ii_3);
                                ds_middle.E_all_linear(:,:,ii)=ds_middle_3.E_all_linear(:,:,ii_3);
                                ds_down.E_all_linear(:,:,ii)=ds_down_3.E_all_linear(:,:,ii_3);
                                ds_middle.alpha_A_vals_linear(:,ii)=ds_middle_3.alpha_A_vals_linear(:,ii_3);
                                ds_down.alpha_A_vals_linear(:,ii)=ds_down_3.alpha_A_vals_linear(:,ii_3);
                                ds_middle.alpha_B_vals_linear(:,ii)=ds_middle_3.alpha_B_vals_linear(:,ii_3);
                                ds_down.alpha_B_vals_linear(:,ii)=ds_down_3.alpha_B_vals_linear(:,ii_3);
                                ds_middle.S_A_vals_linear(:,ii)=ds_middle_3.S_A_vals_linear(:,ii_3);
                                ds_down.S_A_vals_linear(:,ii)=ds_down_3.S_A_vals_linear(:,ii_3);
                                ds_middle.S_B_vals_linear(:,ii)=ds_middle_3.S_B_vals_linear(:,ii_3);
                                ds_down.S_B_vals_linear(:,ii)=ds_down_3.S_B_vals_linear(:,ii_3);
                                ds_middle.Sigma_vals_linear(:,ii)=ds_middle_3.Sigma_vals_linear(:,ii_3);
                                ds_down.Sigma_vals_linear(:,ii)=ds_down_3.Sigma_vals_linear(:,ii_3);
                                ds_middle.E_all_nonlinear(:,:,ii)=ds_middle_3.E_all_nonlinear(:,:,ii_3);
                                ds_down.E_all_nonlinear(:,:,ii)=ds_down_3.E_all_nonlinear(:,:,ii_3);
                                ds_middle.alpha_A_vals_nonlinear(:,ii)=ds_middle_3.alpha_A_vals_nonlinear(:,ii_3);
                                ds_down.alpha_A_vals_nonlinear(:,ii)=ds_down_3.alpha_A_vals_nonlinear(:,ii_3);
                                ds_middle.alpha_B_vals_nonlinear(:,ii)=ds_middle_3.alpha_B_vals_nonlinear(:,ii_3);
                                ds_down.alpha_B_vals_nonlinear(:,ii)=ds_down_3.alpha_B_vals_nonlinear(:,ii_3);
                                ds_middle.S_A_vals_nonlinear(:,ii)=ds_middle_3.S_A_vals_nonlinear(:,ii_3);
                                ds_down.S_A_vals_nonlinear(:,ii)=ds_down_3.S_A_vals_nonlinear(:,ii_3);
                                ds_middle.S_B_vals_nonlinear(:,ii)=ds_middle_3.S_B_vals_nonlinear(:,ii_3);
                                ds_down.S_B_vals_nonlinear(:,ii)=ds_down_3.S_B_vals_nonlinear(:,ii_3);
                                ds_middle.Sigma_vals_nonlinear(:,ii)=ds_middle_3.Sigma_vals_nonlinear(:,ii_3);
                                ds_down.Sigma_vals_nonlinear(:,ii)=ds_down_3.Sigma_vals_nonlinear(:,ii_3);
                                ds_middle.phi_end_nonlinear(:,ii)=ds_middle_3.phi_end_nonlinear(:,ii_3);
                                ds_down.phi_end_nonlinear(:,ii)=ds_down_3.phi_end_nonlinear(:,ii_3);
                            elseif rem(sortedRvec(rr),2)==0&&rem(ss,2)==0
                                ii_4=ii_4+1;
                                ds_middle.parameter_set(ii,:)=ds_middle_4.parameter_set(ii_4,:);
                                ds_down.parameter_set(ii,:)=ds_down_4.parameter_set(ii_4,:);
                                ds_middle.E_all_toroid(:,:,ii)=ds_middle_4.E_all_toroid(:,:,ii_4);
                                ds_down.E_all_toroid(:,:,ii)=ds_down_4.E_all_toroid(:,:,ii_4);
                                ds_middle.alpha_A_vals_toroid(:,ii)=ds_middle_4.alpha_A_vals_toroid(:,ii_4);
                                ds_down.alpha_A_vals_toroid(:,ii)=ds_down_4.alpha_A_vals_toroid(:,ii_4);
                                ds_middle.alpha_B_vals_toroid(:,ii)=ds_middle_4.alpha_B_vals_toroid(:,ii_4);
                                ds_down.alpha_B_vals_toroid(:,ii)=ds_down_4.alpha_B_vals_toroid(:,ii_4);
                                ds_middle.S_A_vals_toroid(:,ii)=ds_middle_4.S_A_vals_toroid(:,ii_4);
                                ds_down.S_A_vals_toroid(:,ii)=ds_down_4.S_A_vals_toroid(:,ii_4);
                                ds_middle.S_B_vals_toroid(:,ii)=ds_middle_4.S_B_vals_toroid(:,ii_4);
                                ds_down.S_B_vals_toroid(:,ii)=ds_down_4.S_B_vals_toroid(:,ii_4);
                                ds_middle.Sigma_vals_toroid(:,ii)=ds_middle_4.Sigma_vals_toroid(:,ii_4);
                                ds_down.Sigma_vals_toroid(:,ii)=ds_down_4.Sigma_vals_toroid(:,ii_4);
                                ds_middle.rho_vals(:,ii)=ds_middle_4.rho_vals(:,ii_4);
                                ds_down.rho_vals(:,ii)=ds_down_4.rho_vals(:,ii_4);
                                ds_middle.E_all_linear(:,:,ii)=ds_middle_4.E_all_linear(:,:,ii_4);
                                ds_down.E_all_linear(:,:,ii)=ds_down_4.E_all_linear(:,:,ii_4);
                                ds_middle.alpha_A_vals_linear(:,ii)=ds_middle_4.alpha_A_vals_linear(:,ii_4);
                                ds_down.alpha_A_vals_linear(:,ii)=ds_down_4.alpha_A_vals_linear(:,ii_4);
                                ds_middle.alpha_B_vals_linear(:,ii)=ds_middle_4.alpha_B_vals_linear(:,ii_4);
                                ds_down.alpha_B_vals_linear(:,ii)=ds_down_4.alpha_B_vals_linear(:,ii_4);
                                ds_middle.S_A_vals_linear(:,ii)=ds_middle_4.S_A_vals_linear(:,ii_4);
                                ds_down.S_A_vals_linear(:,ii)=ds_down_4.S_A_vals_linear(:,ii_4);
                                ds_middle.S_B_vals_linear(:,ii)=ds_middle_4.S_B_vals_linear(:,ii_4);
                                ds_down.S_B_vals_linear(:,ii)=ds_down_4.S_B_vals_linear(:,ii_4);
                                ds_middle.Sigma_vals_linear(:,ii)=ds_middle_4.Sigma_vals_linear(:,ii_4);
                                ds_down.Sigma_vals_linear(:,ii)=ds_down_4.Sigma_vals_linear(:,ii_4);
                                ds_middle.E_all_nonlinear(:,:,ii)=ds_middle_4.E_all_nonlinear(:,:,ii_4);
                                ds_down.E_all_nonlinear(:,:,ii)=ds_down_4.E_all_nonlinear(:,:,ii_4);
                                ds_middle.alpha_A_vals_nonlinear(:,ii)=ds_middle_4.alpha_A_vals_nonlinear(:,ii_4);
                                ds_down.alpha_A_vals_nonlinear(:,ii)=ds_down_4.alpha_A_vals_nonlinear(:,ii_4);
                                ds_middle.alpha_B_vals_nonlinear(:,ii)=ds_middle_4.alpha_B_vals_nonlinear(:,ii_4);
                                ds_down.alpha_B_vals_nonlinear(:,ii)=ds_down_4.alpha_B_vals_nonlinear(:,ii_4);
                                ds_middle.S_A_vals_nonlinear(:,ii)=ds_middle_4.S_A_vals_nonlinear(:,ii_4);
                                ds_down.S_A_vals_nonlinear(:,ii)=ds_down_4.S_A_vals_nonlinear(:,ii_4);
                                ds_middle.S_B_vals_nonlinear(:,ii)=ds_middle_4.S_B_vals_nonlinear(:,ii_4);
                                ds_down.S_B_vals_nonlinear(:,ii)=ds_down_4.S_B_vals_nonlinear(:,ii_4);
                                ds_middle.Sigma_vals_nonlinear(:,ii)=ds_middle_4.Sigma_vals_nonlinear(:,ii_4);
                                ds_down.Sigma_vals_nonlinear(:,ii)=ds_down_4.Sigma_vals_nonlinear(:,ii_4);
                                ds_middle.phi_end_nonlinear(:,ii)=ds_middle_4.phi_end_nonlinear(:,ii_4);
                                ds_down.phi_end_nonlinear(:,ii)=ds_down_4.phi_end_nonlinear(:,ii_4);
                            end
                        end
                    end
                end
            end
        end
    end
end



mu = 1.5e-3; % viscosity, Pa.s
sr = 5000; % shear rate, 1/s

kv = 1;
ev = 1; %2e-5
pv = 1;
av = 1;

%% look for sig
% mu = 3e-3; % viscosity, Pa.s
sv = 1:8;
A = [];

% sig_orig_vals = [0.5, 3, 15]*1e-2;
sig_orig_vals = [0.002, 0.003,0.005];
% sig_orig_vals = [];


tic
for rv = 1:23
    %for rv = 1:68

    slice = get_slice(rv,sv,kv,ev, pv, av);
    outputs_analyse = analyse_phi_curve_updown(ds_middle,ds_down,slice, ...
        'MovmedianRange', 5,'Movmedian_Threshold',0.1,'NSplinePoints',1000,...
        'PlotCurves', false);
    ps = get_param_slices(ds_middle,slice);
    sig_vals = ps.sigma_slice;
    R_vals(rv) = ps.R_slice(1);

    % for each value there, re-run to find hphi
    aB = cell2mat({cell2mat(outputs_analyse.all_global_minima).alpha_B_min});
    aA = cell2mat({cell2mat(outputs_analyse.all_global_minima).alpha_A_min});
    Et = cell2mat({cell2mat(outputs_analyse.all_global_minima).E_total_min})' ...
        - ps.kD_slice/2.*(ps.d_slice.^2.*ps.alpha_i_slice.^2./(ps.alpha_i_slice+1));
    ps.epsilon_slice
    ps.kappa_slice
    for ii=1:length(slice)
        R = ps.R_slice(ii);
        d = ps.d_slice(ii);
        kD = ps.kD_slice(ii);
        ai = ps.alpha_i_slice(ii);
        phi_min(rv,ii) = outputs_analyse.phi_at_min(ii);
        %     E_loc_min = outputs_analyse.E_local_min(ii);
        % %     if ~isempty(phi_loc_min{1})
        % %         phi = phi_loc_min{1};
        % %         phi = phi(1);
        % %         Em = E_loc_min{1};
        % %         E_min(ii) = Em(1);
        % %     else
        % %         hphi(ii) = nan;
        % %         E_min(ii) = nan;
        % %         continue
        % %     end
        %     E_min(ii) = Et(ii);
        %     phi = outputs_analyse.phi_at_min(ii)+1e-5;
        %     kappa = ps.kappa_slice(ii);
        %     Sigma = ps.kD_slice(ii)*aB(ii);
        % %     out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, false);
        %     dxphi(ii) = R*cos(phi);
        %     A(ii) = 2*pi*R^2*(1+cos(phi));


        phiv = ds_down.phi_vals;
        y = ds_down.E_all_nonlinear(1,:,slice(ii)) ...
            - kD/2.*(d.^2.*ai.^2./(ai+1));
        dEdphi = diff(y)./diff(phiv);
        phivm = (phiv(1:end-1)+phiv(2:end))/2;
        xm = 2*R.*cos(phivm);

        %     if ii==12
        %         figure();
        %         hold on
        %         plot(phivm,-dEdphi./xm*1e-6)
        %     %     ylim([-1e-12,1e-12])
        %         xlim([0,pi/3]);
        %         figure();
        %         plot(phiv, y)
        %     end

        if length(-dEdphi(phivm<=outputs_analyse.phi_at_min(ii)))<1 %when static phi is zero
            F_max(ii)=-1;
            A(ii) = 2*pi*R^2*(1+cos(phivm(1)));
            phi_max(rv,ii) = phivm(1);
        else
            [F_max(ii), index] = max((-dEdphi(phivm<=outputs_analyse.phi_at_min(ii))./xm(phivm<=outputs_analyse.phi_at_min(ii))*1e-6)); %J/m
            A(ii) = 2*pi*R^2*(1+cos(phivm(index)));
            phi_max(rv,ii) = phivm(index);
        end

    end
    sr_detach = (F_max./(A*mu*1e-12))';
    sr_detach(sr_detach<0) = nan;
    % check min and max to see if it's bounded
    if max(sr_detach)<sr||all(isnan(sr_detach)) % physiological shear rate is too high for the smallest phi we investigate or system is non-adherent statically
        sig_detach(rv) = 0;
        phi_detach(rv) = 0;
    elseif min(sr_detach)>sr %physiological shear rate is too weak compared to all critical shear rates we examined even at high concentrations
        sig_detach(rv) = max(ps.sigma_slice);
        phi_detach(rv) = min(outputs_analyse.phi_at_min);
    else
        index = find(sr_detach<sr, 1);
        if index==1
            index = 2;
        end
        sig_detach(rv) = 10^interp1(sr_detach(index-1:index), ...
            log10(sig_vals(index-1:index)), sr);
        phi_detach(rv) = 10^interp1(sr_detach(index-1:index),...
            log10(phi_max(rv,index-1:index)), sr);
    end

end
toc


% %% percent detached
figure%('Position',[488,240,689,521]);
hold on
xlabel('Particle Diameter ($\mu$m)');
ylabel('Fraction Remaining');

% experimental data
x_exp = [100,220,450,830]/1000;
y_exp = [0.15, 0.83, 0.86, 0.15];
x = x_exp;
y2 = y_exp;
dy = [0.23, 0.91, 0.89, 0.18]-y2;
%sig_detach
clear err
for ii = 1:length(sig_orig_vals)
    sig_orig = sig_orig_vals(ii);
    y = sig_detach/sig_orig;
    y(y>1) = 1;
    %y3 = y([7,21,30,39]);
    %err(ii) = sqrt(sum((y2-y3).^2));
    %for rv = 1:68
    for rv = 1:23
        phi_orig(ii,rv) = interp1(log10(sig_vals), phi_min(rv,:)', log10(sig_orig));
        phi_new(ii,rv) = interp1(log10(sig_vals), phi_max(rv,:)', log10(sig_detach(rv)));
        if phi_new(ii,rv)>phi_orig(ii,rv) && y(rv)<=1
            y(rv)=nan;
        end
    end


    linewidth = 2;
    %if ii==2
    %    linewidth = 4;
    %end
    plot(R_vals*2, y, '-', 'LineWidth',linewidth, 'Color', newcolors(ii),...
        'DisplayName',sprintf('$\\sigma_0 = %0.3g$', sig_orig));
    %     plot(R_vals*2, y, 'o', 'LineWidth',linewidth, 'Color', newcolors(ii),...
    %         'DisplayName',sprintf('$\\sigma_0 = %0.3g$', sig_orig));

    % get diffs for each sig0
    ydiff = y;
    ydiff(isnan(ydiff)) = 0;
    ynearest = interp1(R_vals*2, ydiff, x_exp, 'nearest');
    diff_sim_exp = sum((ynearest-y_exp).^2);
end

errorbar(x,y2,dy, 'o', 'LineWidth', 2,'MarkerSize',12,...
    'MarkerFaceColor','b', 'Color','b', 'DisplayName','Experiment');
legend
%
% a1 = annotation('textbox', 'String',...
%     [sprintf('$\\epsilon n_0 = %0.3g$ $($pJ$/\\mu$m$^2)$ \n', ps.epsilon_slice(1)),...
%     sprintf('$\\alpha_i = %0.3g$ \n', ps.alpha_i_slice(1)),...
%     sprintf('$\\kappa = %0.3g$ pJ \n', ps.kappa_slice(1)),...
%     sprintf('$k_D = %0.3g$ pJ/$\\mu$m$^2$ \n', ps.kD_slice(1)),...
%     sprintf('$\\mu = %0.3g$ Pa.s \n', mu),...
%     sprintf('$\\dot{\\gamma} = %0.3g$ s$^{-1}$', sr)]);
% a1.Position = [0.7346    0.2016    0.2373    0.1544];
% l1.Position = [0.7257    0.4397    0.2624    0.4453];

%% from sig_0 directly, new truncation method

% sig_orig_vals = [0.5, 3, 15]*1e-2;
sig_orig_vals = [0.025, 0.0075,0.004];
% sig_orig_vals = [0.110.12,0.13];

tic
for ii = 1:length(sig_orig_vals)
    sig0 = sig_orig_vals(ii);
    for rv = 1:23
        %     for rv = 10:20
        count = 0;
        for sv = 1:8
            count = count+1;
            slice = get_slice(rv,sv,kv,ev, pv, av);
            outputs_analyse = analyse_phi_curve_updown(ds_middle,ds_down,...
                slice, 'MovmedianRange', 5,'Movmedian_Threshold',0.1,...
                'NSplinePoints',1000,'PlotCurves', false);

            ps = get_param_slices(ds_middle,slice);
            R_vals(rv) = ps.R_slice(1);
            sig_val(sv) = ps.sigma_slice;
            phiv = ds_down.phi_vals;
            x = ps.R_slice.*cos(phiv);
            y = ds_down.E_all_nonlinear(1,:,slice) ...
                - ps.kD_slice/2.*(ps.d_slice.^2.*...
                ps.alpha_i_slice.^2./(ps.alpha_i_slice+1));
            A = 2*pi*ps.R_slice.^2.*(1+cos(phiv));
            F_shear = (A*mu*1e-12*sr);
            dEdphi = diff(y)./diff(phiv);
            phivm = (phiv(1:end-1)+phiv(2:end))/2;
            xm = 2*ps.R_slice.*cos(phivm);
            F_ad = -dEdphi./xm*1e-6;

            phi1(sv) = outputs_analyse.phi_at_min(1);
            [phi2(sv),val2,~,~] = fminbnd(...
                @(phi) -getAdhesionForce(phi, phivm, F_ad), 0, 1.2);
            val2 = -val2;
            F_shear_at_phi2 = interp1(phiv, F_shear, phi2(sv), 'spline');
            diff_at_phi2(sv) = val2 - F_shear_at_phi2;
            if val2>F_shear_at_phi2
                phi3_cross(sv) = true;
                eqFun = @(phi) getShearForce(phi, ps.R_slice(1), mu, sr)...
                    -getAdhesionForce(phi, phivm, F_ad);
                % find minimum
                meqFun(sv) = fminbnd(eqFun, 0, 1.05);
                valeqFun(sv) = eqFun(meqFun(sv));
                % figure
                %                 figure();
                %                 plot(phiv, eqFun(phiv))
                %                 xlim([0,1.2]);
                %                 ylim([valeqFun(sv)*2, -valeqFun(sv)*2])

                phi4(sv) = fzero(eqFun, [0,meqFun(sv)]);
                if eqFun(1.05)>0
                    phi5(sv) = fzero(eqFun, [meqFun(sv),1.05]);
                else
                    phi5(sv) = 1.05;
                end
            else
                phi3_cross(sv) = false;
                phi4(sv) = 0;
                phi5(sv) = 0;
            end

            %             if count == 1
            %                 f1 = figure();
            %                 hold on
            %                 % plot(phiv, -y./x*1e-6);
            %
            %                 plot(phiv, F_shear, 'DisplayName', '$F_\mathrm{drag}$');
            %                 plot(phivm, F_ad, 'DisplayName',...
            %                     sprintf('$\\sigma = %0.3g$', sig_val));
            %                 xlim([0,1.3])
            %                 xlabel('$\phi$')
            %                 ylabel('$F$')
            %                 a1 = annotation('textbox', 'String',...
            %                 [sprintf('$R = %0.3g$ $\\mu$m ', ps.R_slice(1))]);
            %
            %                 % ylim([min(F_shear), max(F_shear)])
            %                 [val, index] = max(F_ad(1:50));
            %                 ylim([0,val*1.2]);
            %             else
            %                 figure(f1);
            %                 plot(phivm, F_ad, 'DisplayName',...
            %                     sprintf('$\\sigma = %0.3g$', sig_val));
            %                 [valn, index] = max(F_ad(1:50));
            %                 if valn>val
            %                     ylim([0,val*1.2]);
            %                 end
            %             end

        end
        % get phi1 for sig_0 by interpolation, must be shape-preserving, so
        % nearest, linear, pchip etc, otherwise next step doesn't work
        phi1_R(rv) = interp1(sig_val, phi1, sig0, 'pchip');
        % if it's zero, particle isn't attached - just continue
        if interp1(sig_val, phi1, sig0, 'pchip')==0
            sig_detach(rv) = nan;
            continue
        else
            [pp,~] = csaps(sig_val, phi1, 1-1e-6);
            phi1_R(rv) = fnval(pp,sig0);
            %             figure();
            %             hold on
            %             plot(sig_val, fnval(pp, sig_val));
            %             plot(sig_val, phi1);
            %             phi1_R(rv) = interp1(sig_val, phi1, sig0, 'spline');
        end
        % find phi3 by interpolation, get points either side and pos and
        % neg difference between F_shear and F_detach
        cross_index = find(phi3_cross,1);
        phi3_R(rv) = interp1(diff_at_phi2,phi2, 0);
        sig_phi3_cross(rv) = interp1(phi2, sig_val, phi3_R(rv));


        % first check that the crossing phi3 is less than phi1_R
        [phi5_unique, ia, ic] = unique(phi5);
        sig_val_unique = sig_val(ia);

        % testing figure
        %         figure();
        %         hold on
        %         plot(sig_val_unique, phi5_unique-phi1_R(rv))
        %         axes1 = gca;
        %         axes1.XScale = 'log';

        if phi3_R(rv)<phi1_R(rv)
            sig_detach(rv) = interp1(phi5_unique-phi1_R(rv),...
                sig_val_unique, 0, 'spline');
        else
            sig_detach(rv) = nan;
        end
        % find when phi5>phi1_R
        indexDetach = find(phi5>phi1_R(rv),1);
    end
    %     y = sigmin/sig0;
    y = sig_detach/sig0;
    y(y>1) = 1;
    y(sig_phi3_cross>sig0) = 1;
    linewidth = 2;
    if ii==2
        linewidth = 4;
    end
    % plot(R_vals*2, y, ':', 'LineWidth',linewidth, 'Color', newcolors(ii),...
    %     'DisplayName',sprintf('$\\sigma_0 = %0.3g$ stable', sig0));

    % get diffs for each sig0
    ydiff = y;
    ydiff(isnan(ydiff)) = 0;
    ynearest = interp1(R_vals*2, ydiff, x_exp, 'nearest');
    diff_sim_exp_2 = sum((ynearest-y_exp).^2);
end
toc

%% functions

function slice = get_slice(rv, sv, kv, ev, pv, av)
% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = linspace(0.05,1,68)/2;                 % um
R_vals = sort([R_vals(1:5:end) R_vals([13:5:53])]);
sigma_vals = logspace(-4,-0.3, 192);
sigma_vals= sigma_vals(1:25:end);
kD_vals = [0.3];                        % picoJ/um^2
kappa_vals = [3e-8];         % picoJ
epsilon_vals = -8.75e-6;              % picoJ/um^2
n0_vals = 1;                                    % fraction
alpha_i_vals = 1e-6/(kD_vals);                            % fraction

ii = 0;
slice = [];
for rr = 1:length(R_vals)
    for ss = 1:length(sigma_vals)
        for kk = 1:length(kD_vals)
            for ee = 1:length(epsilon_vals)
                for nn = 1:length(n0_vals)
                    for pp = 1:length(kappa_vals)
                        for aa = 1:length(alpha_i_vals)
                            ii = ii+1;
                            if any(rr==rv) && any(ss==sv) && ...
                                    any(kk==kv) && any(ee==ev) &&...
                                    any(rr==rv) && any(pp==pv) && any(aa==av)
                                slice = [slice, ii];
                            end
                        end
                    end
                end
            end
        end
    end
end

end

function F = getShearForce(phi, R, mu, sr)

A = 2*pi*R.^2.*(1+cos(phi));
F = (A*mu*1e-12*sr);

end

function F = getAdhesionForce(phi, phi_vals, Force_vals)
% essentially just fits a spline to the data and interpolates

F = interp1(phi_vals, Force_vals, phi, 'spline');

end

function err = diffAtPhi2(phi, phivm, F_ad, R, mu, sr)

% firstly find phi2
[phi2,F_ad_phi2,~,~] = fminbnd(...
    @(phi) -getAdhesionForce(phi, phivm, F_ad), 0, 1.2);
% now find force at that phi2
F_shear_phi2 = getShearForce(phi2, R, mu, sr);
% now get difference between them
err = F_ad_phi2-F_shear_phi2;

end


function err = getErrorSig0(sig0, ds_middle, ds_down, mu, sr)

kv = 1;
ev = 1;
pv = 1;
av = 1;
tic
countrv = 0;
for rv = [6,7,21,30,31,38,39]
    countrv = countrv+1;
    %     for rv = 20:25
    count = 0;
    for sv = 1:96
        count = count+1;
        slice = get_slice(rv,sv,kv,ev, pv, av);
        outputs_analyse = analyse_phi_curve_updown(ds_middle,ds_down,...
            slice, 'MovmedianRange', 5,'Movmedian_Threshold',0.1,...
            'NSplinePoints',1000,'PlotCurves', false);

        ps = get_param_slices(ds_middle,slice);p
        sig_val(sv) = ps.sigma_slice;
        phiv = ds_down.phi_vals;
        x = ps.R_slice.*cos(phiv);
        y = ds_down.E_all_nonlinear(1,:,slice) ...
            - ps.kD_slice/2.*(ps.d_slice.^2.*...
            ps.alpha_i_slice.^2./(ps.alpha_i_slice+1));
        A = 2*pi*ps.R_slice.^2.*(1+cos(phiv));
        F_shear = (A*mu*1e-12*sr);
        dEdphi = diff(y)./diff(phiv);
        phivm = (phiv(1:end-1)+phiv(2:end))/2;
        xm = 2*ps.R_slice.*cos(phivm);
        F_ad = -dEdphi./xm*1e-6;

        eqFun = @(phi) getShearForce(phi, ps.R_slice(1), mu, sr)...
            -getAdhesionForce(phi, phivm, F_ad);
        [minp(sv),~,~,~] = fzero(eqFun, 0);
        if minp(sv)>1.2||minp(sv)<0
            minp(sv) = nan;
        end

        [maxp(sv),~,~,~] = fminbnd(...
            @(phi) -getAdhesionForce(phi, phivm, F_ad), 0, 1.2);
    end
    % get phimin
    phimax = interp1(sig_val, maxp, sig0);
    indexes = find(minp<phimax);
    index = max(indexes);
    %         indexes = [indexes,max(indexes)+1];
    if ~isempty(indexes)
        if ~isnan(minp(index+1))
            sigmin(countrv) = interp1(minp(index:index+1),...
                sig_val(index:index+1), phimax);
        else
            sigmin(countrv) = interp1([minp(index),maxp(index+1)],...
                sig_val(index:index+1), phimax);
        end
        if isnan(sigmin(countrv))
            sigmin(countrv) = sig_val(max(indexes));
        end
        %             sigmin(rv) = sig_val(max(indexes));
    else
        sigmin(countrv) = nan;
    end

    R_vals(countrv) = ps.R_slice;
end
toc
y = sigmin/sig0;
y(y>1) = 1;

y(isnan(y)) = 0;
x = R_vals*2;

x_exp = [100,220,450,830]/1000;
y_exp = [0.15, 0.83, 0.86, 0.15];

y_sim = interp1(x,y,x_exp);
err = sum((y_sim-y_exp).^2);

end



































































