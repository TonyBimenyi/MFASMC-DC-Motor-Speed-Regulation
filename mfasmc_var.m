clc; clear;

% Parameters
rho = 7.5;            % Increase gain for better control response
eta = 1;              % Increase to make the adaptive update faster
lamda = 350;          % Model-free parameter
mu = 0.005;           % Adaptive parameter
epsilon = 1e-5;       % Small threshold for stability
alpha = 15;           % Smoothing factor for error dynamics
T = 0.1;              % Sampling time (kept for calculations but not used in time vector)
gamma1 = 0.45;        % Adjust control gains for faster tracking
gamma2 = 0.15;        % Adjust control gains
gamma3 = 0.45;        % Adjust control gains
gamma4 = 0.15;        % Adjust control gains
beta = 10;            % Sliding mode coefficient
sigma = 95;           % Sliding mode parameter
tau = 1e-5;           % Small damping term

rT = 1024;            % Sample rate
L = 200;              % Time steps
m = 200;              % Data size
n = 600;              % Data size

% Initialization
phi1 = zeros(m+1,1); 
phi2 = zeros(m+1,1); 
phi3 = zeros(m+1,1); 
phi4 = zeros(m+1,1);
mfa1 = zeros(m+1,1); 
mfa2 = zeros(m+1,1); 
mfa3 = zeros(m+1,1); 
mfa4 = zeros(m+1,1);
sm1 = zeros(m,1);  
sm2 = zeros(m,1);  
sm3 = zeros(m,1);  
sm4 = zeros(m,1);
y1 = zeros(m+1,1); 
y2 = zeros(m+1,1); 
y3 = zeros(m+1,1); 
y4 = zeros(m+1,1);
u1 = zeros(m,1);   
u2 = zeros(m,1);   
u3 = zeros(m,1);   
u4 = zeros(m,1);
xi1 = zeros(m,1);  
xi2 = zeros(m,1);  
xi3 = zeros(m,1);  
xi4 = zeros(m,1);
s1 = zeros(m,1);   
s2 = zeros(m,1);   
s3 = zeros(m,1);   
s4 = zeros(m,1);
yd = zeros(m+1,1);

% Desired signal (Reference trajectory)
for k = 1:1:m+1
    % yd(k) = 0.6; % Time-invariant constant reference signal
    yd(k) = 0.6 * sin(0.05 * pi * k) + 0.6 * cos(0.03 * pi * k);
end

for k = 1:m
    % Adaptive Gain update
    if k == 1
        phi1(k) = 4.0; 
        phi2(k) = 4.0; 
        phi3(k) = 4.0; 
        phi4(k) = 4.0;
    elseif k == 2
        phi1(k) = phi1(k-1) + (eta * u1(k-1) / (mu + u1(k-1)^2)) * (y1(k) - phi1(k-1)*u1(k-1));
        phi2(k) = phi2(k-1) + (eta * u2(k-1) / (mu + u2(k-1)^2)) * (y2(k) - phi2(k-1)*u2(k-1));
        phi3(k) = phi3(k-1) + (eta * u3(k-1) / (mu + u3(k-1)^2)) * (y3(k) - phi3(k-1)*u3(k-1));
        phi4(k) = phi4(k-1) + (eta * u4(k-1) / (mu + u4(k-1)^2)) * (y4(k) - phi4(k-1)*u4(k-1));
    else
        phi1(k) = phi1(k-1) + (eta * (u1(k-1) - u1(k-2)) / (mu + (u1(k-1) - u1(k-2))^2)) * (y1(k) - y1(k-1) - phi1(k-1) * (u1(k-1) - u1(k-2)));
        phi2(k) = phi2(k-1) + (eta * (u2(k-1) - u2(k-2)) / (mu + (u2(k-1) - u2(k-2))^2)) * (y2(k) - y2(k-1) - phi2(k-1) * (u2(k-1) - u2(k-2)));
        phi3(k) = phi3(k-1) + (eta * (u3(k-1) - u3(k-2)) / (mu + (u3(k-1) - u3(k-2))^2)) * (y3(k) - y3(k-1) - phi3(k-1) * (u3(k-1) - u3(k-2)));
        phi4(k) = phi4(k-1) + (eta * (u4(k-1) - u4(k-2)) / (mu + (u4(k-1) - u4(k-2))^2)) * (y4(k) - y4(k-1) - phi4(k-1) * (u4(k-1) - u4(k-2)));
    end

    % Stability protection
    if k > 2 && (abs(phi1(k)) <= epsilon || abs(u1(k-1) - u1(k-2)) <= epsilon || sign(phi1(k)) ~= sign(phi1(1)))
        phi1(k) = phi1(1);
    end
    if k > 2 && (abs(phi2(k)) <= epsilon || abs(u2(k-1) - u2(k-2)) <= epsilon || sign(phi2(k)) ~= sign(phi2(1)))
        phi2(k) = phi2(1);
    end
    if k > 2 && (abs(phi3(k)) <= epsilon || abs(u3(k-1) - u3(k-2)) <= epsilon || sign(phi3(k)) ~= sign(phi3(1)))
        phi3(k) = phi3(1);
    end
    if k > 2 && (abs(phi4(k)) <= epsilon || abs(u4(k-1) - u4(k-2)) <= epsilon || sign(phi4(k)) ~= sign(phi4(1)))
        phi4(k) = phi4(1);
    end

    % Error dynamics
    xi1(k) = yd(k) - 2*y1(k) + y4(k);
    xi2(k) = y1(k) - 2*y2(k) + y3(k);
    xi3(k) = y2(k) + yd(k) - 2*y3(k);
    xi4(k) = y1(k) + y3(k) - 2*y4(k);

    % Fix: Handle k=1 case for sliding surfaces
    if k == 1
        s1(k) = 0;
        s2(k) = 0;
        s3(k) = 0;
        s4(k) = 0;
    else
        s1(k) = alpha * xi1(k) - xi1(k-1);
        s2(k) = alpha * xi2(k) - xi2(k-1);
        s3(k) = alpha * xi3(k) - xi3(k-1);
        s4(k) = alpha * xi4(k) - xi4(k-1);
    end

    % MFAC updates
    if k == 1
        mfa1(k) = 0;
        mfa2(k) = 0;
        mfa3(k) = 0;
        mfa4(k) = 0;
    else
        mfa1(k) = mfa1(k-1) + (rho * phi1(k)) / (lamda + abs(phi1(k)^2)) * xi1(k);
        mfa2(k) = mfa2(k-1) + (rho * phi2(k)) / (lamda + abs(phi2(k)^2)) * xi2(k);
        mfa3(k) = mfa3(k-1) + (rho * phi3(k)) / (lamda + abs(phi3(k)^2)) * xi3(k);
        mfa4(k) = mfa4(k-1) + (rho * phi4(k)) / (lamda + abs(phi4(k)^2)) * xi4(k);
    end

    % SMC updates
    if k == 1 
        sm1(k) = 0;
        sm2(k) = 0;
        sm3(k) = 0;
        sm4(k) = 0;
    else
        sm1(k) = sm1(k-1) + (beta * phi1(k)) / (sigma + (phi1(k))^2) * ...
            ( (xi1(k) + (y4(k) - y4(k-1)) + (yd(k+1) - yd(k))) / (1 + 1) ...
            - (xi1(k)) / (alpha * (2)) + tau * sign(s1(k)) );

        sm2(k) = sm2(k-1) + (beta * phi2(k)) / (sigma + (phi2(k)^2)) * ...
            ( (xi2(k) + (y1(k) - y1(k-1)) + (y3(k) - y3(k-1))) / (1 + 1) ...
            - (xi2(k)) / (alpha * 2) + tau * sign(s2(k)) );

        sm3(k) = sm3(k-1) + (beta * phi3(k)) / (sigma + (phi3(k)^2)) * ...
            ( (xi3(k) + (y2(k) - y2(k-1)) + (yd(k+1) - yd(k))) / (1 + 1) ...
            - (xi3(k)) / (alpha * (2)) + tau * sign(s3(k)) );

        sm4(k) = sm4(k-1) + (beta * phi4(k)) / (sigma + (phi4(k)^2)) * ...
            ( (xi4(k) + (y1(k) - y1(k-1)) + (y3(k) - y3(k-1))) / (1 + 1) ...
            - (xi4(k)) / (alpha * (2)) + tau * sign(s4(k)) );
    end

    % Control signal
    if k == 1
        u1(k) = 0.01;
        u2(k) = 0.01;
        u3(k) = 0.01;
        u4(k) = 0.01;
    else
        u1(k) = mfa1(k) + gamma1 * sm1(k);
        u2(k) = mfa2(k) + gamma2 * sm2(k);
        u3(k) = mfa3(k) + gamma3 * sm3(k);
        u4(k) = mfa4(k) + gamma4 * sm4(k);
    end
    if k == 1
        y1(k) = 0.6;
        y2(k) = 0.6;
        y3(k) = 0.6;
        y4(k) = 0.6;
    end

    % Plant model update with nonlinear term and feedforward
    a = 0.5;
    b1 = 1.2 * n / (rT * 0.2);
    b2 = 1.2 * n / (rT * 0.2);
    b3 = 1.2 * n / (rT * 0.2);
    b4 = 1.2 * n / (rT * 0.2);
    
    nonlinearity1 = 0.03; % Coefficient for cubic nonlinearity
    nonlinearity2 = 0.01; % Coefficient for cubic nonlinearity
    nonlinearity3 = 0.02; % Coefficient for cubic nonlinearity
    nonlinearity4 = 0.01; % Coefficient for cubic nonlinearity
    ff_gain = 0.45; % Feedforward gain
    
    % Add cubic nonlinearity and feedforward term
    y1(k+1) = a * y1(k) + b1 * u1(k) - nonlinearity1 * y1(k)^3 + ff_gain;
    y2(k+1) = a * y2(k) + b2 * u2(k) - nonlinearity2 * y2(k)^2 + ff_gain;
    y3(k+1) = a * y3(k) + b3 * u3(k) - nonlinearity3 * y3(k)^3 + ff_gain;
    y4(k+1) = a * y4(k) + b4 * u4(k) - nonlinearity4 * y4(k)^2 + ff_gain;
end

% Plotting
% Time vector for plotting
t = 1:1:m+1;  % Gives 201 points, as expected
font_size = 14;

% Verify time vector and array lengths
disp(['Length of t: ', num2str(length(t))]); % Should display 201
disp(['Length of y1: ', num2str(length(y1))]); % Should be 201
disp(['Length of yd: ', num2str(length(yd))]); % Should be 201
disp(['Length of xi1: ', num2str(length(xi1))]); % Should be 200

% Plot 4 subplots
figure('Position', [100, 100, 15*100, 7.5*100]); % [left, bottom, width, height] in pixels
subplot(2,2,1);
plot(t, yd, '--b', 'LineWidth', 2.5); hold on;
plot(t, y1, '-.g', 'LineWidth', 2.5);
title('Agent 1'); grid off;
legend('y_d(k)','y_1(k)','Orientation', 'horizontal');
set(gca, 'FontSize', font_size);
xlim([0 m]); % X-axis starts from 0
ylim([-1.5 3.5]); % Y-axis limits for Agent 1

zoom_x_start = 80; % Start of zoomed x-range
zoom_x_end = 100; % End of zoomed x-range
axes('Position', [0.20,0.79,0.13,0.10]);
box on; hold on;
plot(t, yd, '--b', 'LineWidth', 2.5);
plot(t, y1, '-.g', 'LineWidth', 2.5);
xlim([zoom_x_start zoom_x_end]);
yticks([-0.4,0,0.2]);
set(gca, 'FontSize', font_size);

subplot(2,2,2);
plot(t, yd, '--b', 'LineWidth', 2.5); hold on;
plot(t, y2, '-.g', 'LineWidth', 2.5);
title('Agent 2'); grid off;
legend('y_d(k)','y_2(k)','Orientation', 'horizontal');
set(gca, 'FontSize', font_size);
xlim([0 m]); % X-axis starts from 0
ylim([-1.5 3.5]); % Y-axis limits for Agent 2

axes('Position', [0.64,0.79,0.13,0.10]);
box on; hold on;
plot(t, yd, '--b', 'LineWidth', 2.5);
plot(t, y2, '-.g', 'LineWidth', 2.5);
xlim([zoom_x_start zoom_x_end]);
yticks([-0.4,0,0.2]);
set(gca, 'FontSize', font_size);

subplot(2,2,3);
plot(t, yd, '--b', 'LineWidth', 2.5); hold on;
plot(t, y3, '-.g', 'LineWidth', 2.5);
title('Agent 3'); grid off;
legend('y_d(k)','y_3(k)','Orientation', 'horizontal');
set(gca, 'FontSize', font_size);
xlim([0 m]); % X-axis starts from 0
ylim([-1.5 3.5]); % Y-axis limits for Agent 3

axes('Position', [0.20,0.315,0.13,0.10]);
box on; hold on;
plot(t, yd, '--b', 'LineWidth', 2.5);
plot(t, y3, '-.g', 'LineWidth', 2.5);
xlim([zoom_x_start zoom_x_end]);
yticks([-0.4,0,0.2]);
set(gca, 'FontSize', font_size);

subplot(2,2,4);
plot(t, yd, '--b', 'LineWidth', 2.5); hold on;
plot(t, y4, '-.g', 'LineWidth', 2.5);
title('Agent 4'); grid off;
legend('y_d(k)','y_4(k)','Orientation', 'horizontal');
set(gca, 'FontSize', font_size);
xlim([0 m]); % X-axis starts from 0
ylim([-1.5 3.5]); % Y-axis limits for Agent 4

axes('Position', [0.64,0.315,0.13,0.10]);
box on; hold on;
plot(t, yd, '--b', 'LineWidth', 2.5);
plot(t, y4, '-.g', 'LineWidth', 2.5);
xlim([zoom_x_start zoom_x_end]);
yticks([-0.4,0,0.2]);
set(gca, 'FontSize', font_size);







% Plot 4 subplots


figure('Position', [100, 100, 15*100, 7.5*100]); % [left, bottom, width, height] in pixels
subplot(2,2,1);
plot(t(1:end-1), xi1, '-.g', 'LineWidth', 2.5);
title('Agent 1'); grid off;
legend('\xi_1(k)','y_1','Orientation', 'horizontal');
set(gca, 'FontSize', font_size);
xlim([0 m]); % X-axis starts from 0
% ylim([0.5 0.8]); % Y-axis limits for Agent 1

zoom_x_start_xi = 50; % Start of zoomed x-range
zoom_x_end_xi = 80; % End of zoomed x-range
axes('Position', [0.20,0.75,0.15,0.13]);
box on; hold on;
plot(t(1:end-1), xi1, '-.g', 'LineWidth', 2.5);
xlim([zoom_x_start_xi zoom_x_end_xi]);
% yticks([-2.000000000000000e-04,0,2.000000000000000e-04,4.000000000000001e-04,6.000000000000001e-04]);
set(gca, 'FontSize', font_size);

subplot(2,2,2);
plot(t(1:end-1), xi2, '-.g', 'LineWidth', 2.5);
title('Agent 2'); grid off;
legend('\xi_2(k)','y_2','Orientation', 'horizontal');
set(gca, 'FontSize', font_size);
xlim([0 m]); % X-axis starts from 0
% ylim([0.41 0.99]); % Y-axis limits for Agent 2

axes('Position', [0.65,0.64,0.15,0.13]);
box on; hold on;
plot(t(1:end-1), xi2, '-.g', 'LineWidth', 2.5);
xlim([zoom_x_start_xi zoom_x_end_xi]);
% yticks([-2.000000000000000e-04,0,2.000000000000000e-04,4.000000000000001e-04,6.000000000000001e-04]);
set(gca, 'FontSize', font_size);

subplot(2,2,3);
plot(t(1:end-1), xi2, '-.g', 'LineWidth', 2.5);
title('Agent 3'); grid off;
legend('\xi_3(k)','y_3','Orientation', 'horizontal');
set(gca, 'FontSize', font_size);
xlim([0 m]); % X-axis starts from 0
% ylim([0.5 0.8]); % Y-axis limits for Agent 1

axes('Position', [0.20,0.170,0.15,0.13]);
box on; hold on;
plot(t(1:end-1), xi2, '-.g', 'LineWidth', 2.5);
xlim([zoom_x_start_xi zoom_x_end_xi]);
% yticks([0.599,0.6,0.601]);
set(gca, 'FontSize', font_size);

subplot(2,2,4);
plot(t(1:end-1), xi4, '-.g', 'LineWidth', 2.5);
title('Agent 4'); grid off;
legend('\xi_4(k)','y_4','Orientation', 'horizontal');
set(gca, 'FontSize', font_size);
xlim([0 m]); % X-axis starts from 0
% ylim([0.41 0.99]); % Y-axis limits for Agent 2

axes('Position', [0.65,0.170,0.15,0.13]);
box on; hold on;
plot(t(1:end-1), xi4, '-.g', 'LineWidth', 2.5);
xlim([zoom_x_start_xi zoom_x_end_xi]);
% xticks([64,64.5,65])
% yticks([-1.000000000000000e-04,0,1.000000000000000e-04]);
set(gca, 'FontSize', font_size);
hold off;