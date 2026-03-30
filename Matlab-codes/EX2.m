clc;
clear;
close all;
% Parameters
c = 3e8; % speed of light
fc = 6e9; % frequency
lambda = c / fc; % wavelength
NZ = 16; % number of sensors along z
NY = 16; % number of sensors along y
N_tot = NZ * NY;
d_z = lambda / 2; % sensor spacing along z [m]
d_y = lambda / 2; % sensor spacing along y [m]
theta_1deg = 90; % UE el. angle from 0 to 180
theta_1 = deg2rad(theta_1deg);
phi_1deg = 0; % UE az. from -90 to 90
phi_1 = deg2rad(phi_1deg);
sigma_n2 = 1e-5; % noise variance
% Define the directivity function
directivity_func = @(theta, phi) 0.25 * (1 - cos(2 * theta)) .* (1 + cos(phi));
% Compute the beamforming filter using the kron function
a_theta = @(theta, phi) kron(exp(1j * 2 * pi * d_y / lambda * (0:NY-1)' * sin(phi)), ...
                             exp(1j * 2 * pi * d_z / lambda * (0:NZ-1)' * cos(theta)));
% Generate angle grids
angles_1 = 361;
theta = linspace(0, pi, angles_1);
angles_2 = 721;
phi = linspace(-pi, pi, angles_2);
[PHI, THETA] = meshgrid(phi, theta);
% Generate the UPA steering vectors for each angle pair
A_theta = zeros(N_tot, angles_1, angles_2);
for i = 1:angles_1
    for j = 1:angles_2
        A_theta(:, i, j) = a_theta(THETA(i, j), PHI(i, j));
    end
end
% Compute the UPA pattern for isotropic antenna elements
w = a_theta(theta_1, phi_1); % Beamforming vector for the UE
AF = zeros(angles_1, angles_2);
for i = 1:angles_1
    for j = 1:angles_2
        AF(i, j) = abs(w' * squeeze(A_theta(:, i, j)));
    end
end
% Plot the 3-D UPA pattern with isotropic antenna elements
figure;
mesh(PHI, THETA, AF);
axis equal;
xlabel('Phi');
ylabel('Theta');
zlabel('AF');
title('3-D UPA Pattern with Isotropic Antenna Elements');
% Apply pattern multiplication with the directivity function
D = directivity_func(THETA, PHI);
AF_directive = AF .* D;
% Plot the 3-D UPA pattern with directive antenna elements
figure;
mesh(PHI, THETA, AF_directive);
axis equal;
xlabel('Phi');
ylabel('Theta');
zlabel('AF');
title('3-D UPA Pattern with Directive Antenna Elements');
% Compare the differences
disp('Differences between isotropic and directive antenna elements can be observed in the plots.');
% Test with different angle pairs
angle_pairs = [deg2rad(105), deg2rad(30); deg2rad(70), deg2rad(-45)];
for k = 1:size(angle_pairs, 1)
    theta_test = angle_pairs(k, 1);
    phi_test = angle_pairs(k, 2);
    
    % Beamforming vector
    w_test = a_theta(theta_test, phi_test);
    
    % UPA pattern for isotropic elements
    AF_test = zeros(angles_1, angles_2);
    for i = 1:angles_1
        for j = 1:angles_2
            AF_test(i, j) = abs(w_test' * squeeze(A_theta(:, i, j)));
        end
    end
    
    % Plot for isotropic elements
    figure;
    mesh(PHI, THETA, AF_test);
    axis equal;
    xlabel('Phi');
    ylabel('Theta');
    zlabel('AF');
    title(sprintf('UPA Pattern with Isotropic Antenna Elements (\\theta=%.1f°, \\phi=%.1f°)', rad2deg(theta_test), rad2deg(phi_test)));
    
    % UPA pattern for directive elements
    AF_test_directive = AF_test .* directivity_func(THETA, PHI);
    
    % Plot for directive elements
    figure;
    mesh(PHI, THETA, AF_test_directive);
    axis equal;
    xlabel('Phi');
    ylabel('Theta');
    zlabel('AF');
    title(sprintf('UPA Pattern with Directive Antenna Elements (\\theta=%.1f°, \\phi=%.1f°)', rad2deg(theta_test), rad2deg(phi_test)));
end
% Set NZ = 4 and NY = 32, plot and compare
NZ = 4;
NY = 32;
N_tot = NZ * NY;
d_z = lambda / 2;
d_y = lambda / 2;
a_theta = @(theta, phi) kron(exp(1j * 2 * pi * d_y / lambda * (0:NY-1)' * sin(phi)), ...
                             exp(1j * 2 * pi * d_z / lambda * (0:NZ-1)' * cos(theta)));
A_theta = zeros(N_tot, angles_1, angles_2);
for i = 1:angles_1
    for j = 1:angles_2
        A_theta(:, i, j) = a_theta(THETA(i, j), PHI(i, j));
    end
end
w = a_theta(deg2rad(100), deg2rad(60));
AF = zeros(angles_1, angles_2);
for i = 1:angles_1
    for j = 1:angles_2
        AF(i, j) = abs(w' * squeeze(A_theta(:, i, j)));
    end
end
figure;
mesh(PHI, THETA, AF);
axis equal;
xlabel('Phi');
ylabel('Theta');
zlabel('AF');
title('UPA Pattern with Isotropic Antenna Elements (NZ=4, NY=32, θ=100°, φ=60°)');
AF_directive = AF .* directivity_func(THETA, PHI);
figure;
mesh(PHI, THETA, AF_directive);
axis equal;
xlabel('Phi');
ylabel('Theta');
zlabel('AF');
title('UPA Pattern with Directive Antenna Elements (NZ=4, NY=32, θ=100°, φ=60°)');
% Interferers' DoAs
theta_intfdeg = [86, 85, 80, 100, 105];
phi_intfdeg = [4, 20, 5, -15, 15];
theta_intf = deg2rad(theta_intfdeg);
phi_intf = deg2rad(phi_intfdeg);
% Compute the covariance matrix of interferers
R_i = zeros(N_tot, N_tot);
for k = 1:length(theta_intf)
    a_intf = a_theta(theta_intf(k), phi_intf(k));
    R_i = R_i + (a_intf * a_intf');
end
R_i = R_i + sigma_n2 * eye(N_tot);
% Compute the MVDR beamformer
a_ue = a_theta(theta_1, phi_1);
R_i_inv = inv(R_i);
w_mvdr = (R_i_inv * a_ue) / (a_ue' * R_i_inv * a_ue);
% Compute the UPA pattern for MVDR
AF_mvdr = zeros(angles_1, angles_2);
for i = 1:angles_1
    for j = 1:angles_2
        AF_mvdr(i, j) = abs(w_mvdr' * squeeze(A_theta(:, i, j)));
    end
end
% Apply pattern multiplication with the directivity function
AF_mvdr_directive = AF_mvdr .* directivity_func(THETA, PHI);
% Plot the pattern with isotropic antenna elements (directive=0)
figure;
mesh(PHI, THETA, AF_mvdr);
axis equal;
xlabel('Phi');
ylabel('Theta');
zlabel('AF');
title('MVDR Beamforming Pattern with Isotropic Antenna Elements');
% Plot the pattern with directive antenna elements (directive=1)
figure;
mesh(PHI, THETA, AF_mvdr_directive);
axis equal;
xlabel('Phi');
ylabel('Theta');
zlabel('AF');
title('MVDR Beamforming Pattern with Directive Antenna Elements');
% Change the DoA of the first interferer and re-compute
theta_intf(1) = deg2rad(88);
phi_intf(1) = deg2rad(2);
R_i = zeros(N_tot, N_tot);
for k = 1:length(theta_intf)
    a_intf = a_theta(theta_intf(k), phi_intf(k));
    R_i = R_i + (a_intf * a_intf');
end
R_i = R_i + sigma_n2 * eye(N_tot);
R_i_inv = inv(R_i);
w_mvdr = (R_i_inv * a_ue) / (a_ue' * R_i_inv * a_ue);
AF_mvdr = zeros(angles_1, angles_2);
for i = 1:angles_1
    for j = 1:angles_2
        AF_mvdr(i, j) = abs(w_mvdr' * squeeze(A_theta(:, i, j)));
    end
end
AF_mvdr_directive = AF_mvdr .* directivity_func(THETA, PHI);
% Plot the pattern with the new interferer DoA
figure;
mesh(PHI, THETA, AF_mvdr_directive);
axis equal;
xlabel('Phi');
ylabel('Theta');
zlabel('AF');
title('MVDR Beamforming Pattern with New Interferer DoA (θ=88°, φ=2°)');
% Increase noise variance and analyze the effect
sigma_n2 = 1e3;
R_i = zeros(N_tot, N_tot);
for k = 1:length(theta_intf)
    a_intf = a_theta(theta_intf(k), phi_intf(k));
    R_i = R_i + (a_intf * a_intf');
end
R_i = R_i + sigma_n2 * eye(N_tot);
R_i_inv = inv(R_i);
w_mvdr = (R_i_inv * a_ue) / (a_ue' * R_i_inv * a_ue);
AF_mvdr = zeros(angles_1, angles_2);
for i = 1:angles_1
    for j = 1:angles_2
        AF_mvdr(i, j) = abs(w_mvdr' * squeeze(A_theta(:, i, j)));
    end
end
AF_mvdr_directive = AF_mvdr .* directivity_func(THETA, PHI);
% Plot the pattern with increased noise variance
figure;
mesh(PHI, THETA, AF_mvdr_directive);
axis equal;
xlabel('Phi');
ylabel('Theta');
zlabel('AF');
title('MVDR Beamforming Pattern with Increased Noise Variance (σ_n^2=10^3)');
disp('The increased noise variance affects the pattern by making it less sharp and increasing the sidelobes.');
