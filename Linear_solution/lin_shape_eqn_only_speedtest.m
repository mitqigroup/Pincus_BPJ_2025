N = 1e4;
d = 10;
R = 1;
epsilon = -1;
n0 = 1;
Sigma = -1;
kappa = 2e-5;
lambda = sqrt(kappa/Sigma);

% get the shape of the free region
% phi_vals = linspace(0.001, pi/8, 5);
% h0_vals = linspace(-1.5,0.5,5);
phi_vals = 0.001;
h0_vals = -1.5;

E = zeros(length(phi_vals), length(h0_vals));
E2 = zeros(size(E));
A = zeros(size(E));
A2 = zeros(size(E2));
C = zeros(length(phi_vals), length(h0_vals),4);
C2 = zeros(length(phi_vals), length(h0_vals),4);

tic
for ii = 1:length(phi_vals)
for jj = 1:length(h0_vals)
phi = phi_vals(ii);
h0 = h0_vals(jj);
r_phi = sin(phi)*R;
r = linspace(r_phi, d/2,N);
S_A = 2*pi*R^2*(1-cos(phi));

[C2(ii,jj,:),A2(ii,jj),E2(ii,jj)] = free_shape_linear_no_curve(r_phi, d, phi, kappa, Sigma, h0);

end
end
toc

tic
for ii = 1:length(phi_vals)
for jj = 1:length(h0_vals)

phi = phi_vals(ii);
h0 = h0_vals(jj);
r_phi = sin(phi)*R;
r = linspace(r_phi, d/2,N);
S_A = 2*pi*R^2*(1-cos(phi));

[h,C(ii,jj,:),A(ii,jj),E(ii,jj),hderiv,lap_h] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h0);

end
end
toc

% figure();
% hold on
% plot(r, lap_h.^2)
% 
% figure();
% hold on
% plot(r, cumtrapz(r,lap_h.^2))
% 
% figure();
% hold on
% plot(r, r.*sqrt(1+hderiv.^2))
% 
% figure();
% hold on
% plot(r, cumtrapz(r, r.*sqrt(1+hderiv.^2)))
