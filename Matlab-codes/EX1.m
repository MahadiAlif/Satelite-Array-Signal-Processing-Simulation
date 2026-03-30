clc; clear all; clear memory;

% Step 1: Initialization
fc = 900e6;
c = 3e8;
lambda = c / fc;
theta1_deg = 0;
theta1_rad = deg2rad(theta1_deg);
S = 3601;
theta = linspace(-pi/2, pi/2, S);
m = (0:50).';

% Step 2: Compute beamforming vector for theta1
M = 8;
d = lambda / 2;
m_vec = (0:M-1).';
psi1 = 2 * pi * d / lambda * sin(theta1_rad);
w = exp(1j * m_vec * psi1) / sqrt(M);

% Step 3: Compute steering matrix V and array pattern
V = zeros(M, S);
for i = 1:S
    psi = 2 * pi * d / lambda * sin(theta(i));
    V(:, i) = exp(1j * m_vec * psi);
end
pattern = w' * V;

% Step 4: Plot array factor and diric function
figure;
plot(rad2deg(theta), abs(pattern), 'b', 'LineWidth', 1.5); hold on;
xlabel('\theta (degrees)');
ylabel('Array Factor Magnitude');
title('Beamforming Pattern');
grid on;
psi = 2 * pi * d / lambda * (sin(theta) - sin(theta1_rad));
plot(rad2deg(theta), abs(diric(psi, M)), '--r', 'LineWidth', 1.2);
legend('Beamformer Output', 'Dirichlet Function');

% Polar plot
figure;
polarplot(theta, abs(pattern), 'r', 'LineWidth', 2);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaLim = [-90 90];
title(['Conventional beamforming with M=', num2str(M), ', \theta_1 = ', num2str(theta1_deg), '^o'], 'FontSize', 14);

% Step 5 - Case 1: Vary M
figure;
M_values = [2, 4, 8, 50];
for i = 1:4
    M = M_values(i);
    m = (0:M-1).';
    d = lambda / 2;
    psi1 = 2 * pi * d / lambda * sin(theta1_rad);
    w = exp(1j * m * psi1) / sqrt(M);
    V = zeros(M, S);
    for k = 1:S
        psi = 2 * pi * d / lambda * sin(theta(k));
        V(:, k) = exp(1j * m * psi);
    end
    pattern = w' * V;
    subplot(2, 2, i);
    polarplot(theta, abs(pattern), 'LineWidth', 2);
    ax = gca; ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'clockwise'; ax.ThetaLim = [-90 90];
    title(['M = ', num2str(M), ', d = \lambda/2']);
end

% Step 5 - Case 2: Vary d/lambda
figure;
d_ratios = [0.25, 0.5, 1, 2];
M = 8;
m = (0:M-1).';
for i = 1:4
    d = d_ratios(i) * lambda;
    psi1 = 2 * pi * d / lambda * sin(theta1_rad);
    w = exp(1j * m * psi1) / sqrt(M);
    V = zeros(M, S);
    for k = 1:S
        psi = 2 * pi * d / lambda * sin(theta(k));
        V(:, k) = exp(1j * m * psi);
    end
    pattern = w' * V;
    subplot(2, 2, i);
    polarplot(theta, abs(pattern), 'LineWidth', 2);
    ax = gca; ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'clockwise'; ax.ThetaLim = [-90 90];
    title(['M = 8, d/\lambda = ', num2str(d_ratios(i))]);
end

% Step 8–14: MVDR Beamforming
% Interferers
theta1_deg = 0;
theta1_rad = deg2rad(theta1_deg);
theta_int_deg = [20, -40, 60, -75, 80];
theta_all_deg = [theta1_deg, theta_int_deg];
theta_all_rad = deg2rad(theta_all_deg);
M = 8; d = lambda / 2; m = (0:M-1).';

a = @(th) exp(1j * 2 * pi * d / lambda * m * sin(th));

A = zeros(M, length(theta_all_rad));
for k = 1:length(theta_all_rad)
    A(:,k) = a(theta_all_rad(k));
end

sigma_n = 1e-5;
Ry = A * A' + sigma_n * eye(M);

a1 = a(theta1_rad);
w_mvdr = (Ry \ a1);
w_mvdr = w_mvdr / norm(w_mvdr);

V = zeros(M, S);
for i = 1:S
    V(:, i) = a(theta(i));
end
pattern = w_mvdr' * V;

figure;
polarplot(theta, abs(pattern), 'b', 'LineWidth', 2);
ax = gca; ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'clockwise'; ax.ThetaLim = [-90 90];
title('MVDR Beamforming Pattern with Interferers');
hold on;
rho_max = max(abs(pattern));
for ang = theta_int_deg
    polarplot([deg2rad(ang), deg2rad(ang)], [0, rho_max], '--r');
end
