clc; clear;

% Common Parameters
rho = 7.5;            % Gain for control response
eta = 1;              % Adaptive update speed
lamda = 350;          % Model-free parameter
mu = 0.005;           % Adaptive parameter
epsilon = 1e-5;       % Stability threshold
alpha = 15;           % Smoothing factor for error dynamics
T = 0.1;              % Sampling time
m = 200;              % Time steps
rT = 1024;            % Sample rate
font_size = 14;

% Method 1 (MFAC + SMC) Parameters
gamma1_m1 = 0.45;     % Control gains
gamma2_m1 = 0.25;
gamma3_m1 = 0.45;
gamma4_m1 = 0.25;
beta_m1 = 10;         % Sliding mode coefficient
sigma_m1 = 95;        % Sliding mode parameter
tau_m1 = 1e-5;        % Damping term
n_m1 = 600;           % Data size for Method 1
a_m1 = 0.5;           % Plant parameter
ff_gain_m1 = 0.45;    % Feedforward gain

% Method 2 (MFAC only) Parameters
gamma1_m2 = 0;        % Disable SMC
gamma2_m2 = 0;
gamma3_m2 = 0;
gamma4_m2 = 0;
n_m2 = 1024;          % Data size for Method 2
a_m2 = 0.71;           % Plant parameter
ff_gain_m2 = 0.35;     % Feedforward gain

% Nonlinearity coefficients (same for both methods)
nonlinearity1 = 0.02;
nonlinearity2 = 0.02;
nonlinearity3 = 0.02;
nonlinearity4 = 0.02;

% Initialization for both methods
phi1_m1 = zeros(m+1,1); phi2_m1 = zeros(m+1,1); phi3_m1 = zeros(m+1,1); phi4_m1 = zeros(m+1,1);
mfa1_m1 = zeros(m+1,1); mfa2_m1 = zeros(m+1,1); mfa3_m1 = zeros(m+1,1); mfa4_m1 = zeros(m+1,1);
sm1_m1 = zeros(m,1); sm2_m1 = zeros(m,1); sm3_m1 = zeros(m,1); sm4_m1 = zeros(m,1);
y1_m1 = zeros(m+1,1); y2_m1 = zeros(m+1,1); y3_m1 = zeros(m+1,1); y4_m1 = zeros(m+1,1);
u1_m1 = zeros(m,1); u2_m1 = zeros(m,1); u3_m1 = zeros(m,1); u4_m1 = zeros(m,1);
xi1_m1 = zeros(m,1); xi2_m1 = zeros(m,1); xi3_m1 = zeros(m,1); xi4_m1 = zeros(m,1);
s1_m1 = zeros(m,1); s2_m1 = zeros(m,1); s3_m1 = zeros(m,1); s4_m1 = zeros(m,1);

phi1_m2 = zeros(m+1,1); phi2_m2 = zeros(m+1,1); phi3_m2 = zeros(m+1,1); phi4_m2 = zeros(m+1,1);
mfa1_m2 = zeros(m+1,1); mfa2_m2 = zeros(m+1,1); mfa3_m2 = zeros(m+1,1); mfa4_m2 = zeros(m+1,1);
y1_m2 = zeros(m+1,1); y2_m2 = zeros(m+1,1); y3_m2 = zeros(m+1,1); y4_m2 = zeros(m+1,1);
u1_m2 = zeros(m,1); u2_m2 = zeros(m,1); u3_m2 = zeros(m,1); u4_m2 = zeros(m,1);
xi1_m2 = zeros(m,1); xi2_m2 = zeros(m,1); xi3_m2 = zeros(m,1); xi4_m2 = zeros(m,1);

yd = zeros(m+1,1);

% Desired signal (Reference trajectory)
for k = 1:m+1
    % yd(k) = 0.6 * sin(0.05 * pi * k) + 0.6 * cos(0.03 * pi * k);
    yd(k) = 0.6; % Time-invariant constant reference signalx
end

% Method 1 (MFAC + SMC) Simulation
for k = 1:m
    % Adaptive Gain update
    if k == 1
        phi1_m1(k) = 1; phi2_m1(k) = 1; phi3_m1(k) = 1; phi4_m1(k) = 1;
    elseif k == 2
        phi1_m1(k) = phi1_m1(k-1) + (eta * u1_m1(k-1) / (mu + u1_m1(k-1)^2)) * (y1_m1(k) - phi1_m1(k-1)*u1_m1(k-1));
        phi2_m1(k) = phi2_m1(k-1) + (eta * u2_m1(k-1) / (mu + u2_m1(k-1)^2)) * (y2_m1(k) - phi2_m1(k-1)*u2_m1(k-1));
        phi3_m1(k) = phi3_m1(k-1) + (eta * u3_m1(k-1) / (mu + u3_m1(k-1)^2)) * (y3_m1(k) - phi3_m1(k-1)*u3_m1(k-1));
        phi4_m1(k) = phi4_m1(k-1) + (eta * u4_m1(k-1) / (mu + u4_m1(k-1)^2)) * (y4_m1(k) - phi4_m1(k-1)*u4_m1(k-1));
    else
        phi1_m1(k) = phi1_m1(k-1) + (eta * (u1_m1(k-1) - u1_m1(k-2)) / (mu + (u1_m1(k-1) - u1_m1(k-2))^2)) * (y1_m1(k) - y1_m1(k-1) - phi1_m1(k-1) * (u1_m1(k-1) - u1_m1(k-2)));
        phi2_m1(k) = phi2_m1(k-1) + (eta * (u2_m1(k-1) - u2_m1(k-2)) / (mu + (u2_m1(k-1) - u2_m1(k-2))^2)) * (y2_m1(k) - y2_m1(k-1) - phi2_m1(k-1) * (u2_m1(k-1) - u2_m1(k-2)));
        phi3_m1(k) = phi3_m1(k-1) + (eta * (u3_m1(k-1) - u3_m1(k-2)) / (mu + (u3_m1(k-1) - u3_m1(k-2))^2)) * (y3_m1(k) - y3_m1(k-1) - phi3_m1(k-1) * (u3_m1(k-1) - u3_m1(k-2)));
        phi4_m1(k) = phi4_m1(k-1) + (eta * (u4_m1(k-1) - u4_m1(k-2)) / (mu + (u4_m1(k-1) - u4_m1(k-2))^2)) * (y4_m1(k) - y4_m1(k-1) - phi4_m1(k-1) * (u4_m1(k-1) - u4_m1(k-2)));
    end

    % Stability protection
    if k > 2 && (abs(phi1_m1(k)) <= epsilon || abs(u1_m1(k-1) - u1_m1(k-2)) <= epsilon || sign(phi1_m1(k)) ~= sign(phi1_m1(1)))
        phi1_m1(k) = phi1_m1(1);
    end
    if k > 2 && (abs(phi2_m1(k)) <= epsilon || abs(u2_m1(k-1) - u2_m1(k-2)) <= epsilon || sign(phi2_m1(k)) ~= sign(phi2_m1(1)))
        phi2_m1(k) = phi2_m1(1);
    end
    if k > 2 && (abs(phi3_m1(k)) <= epsilon || abs(u3_m1(k-1) - u3_m1(k-2)) <= epsilon || sign(phi3_m1(k)) ~= sign(phi3_m1(1)))
        phi3_m1(k) = phi3_m1(1);
    end
    if k > 2 && (abs(phi4_m1(k)) <= epsilon || abs(u4_m1(k-1) - u4_m1(k-2)) <= epsilon || sign(phi4_m1(k)) ~= sign(phi4_m1(1)))
        phi4_m1(k) = phi4_m1(1);
    end

    % Error dynamics
    xi1_m1(k) = yd(k) - 2*y1_m1(k) + y4_m1(k);
    xi2_m1(k) = y1_m1(k) - 2*y2_m1(k) + y3_m1(k);
    xi3_m1(k) = y2_m1(k) + yd(k) - 2*y3_m1(k);
    xi4_m1(k) = y1_m1(k) + y3_m1(k) - 2*y4_m1(k);

    % Sliding surfaces
    if k == 1
        s1_m1(k) = 0; s2_m1(k) = 0; s3_m1(k) = 0; s4_m1(k) = 0;
    else
        s1_m1(k) = alpha * xi1_m1(k) - xi1_m1(k-1);
        s2_m1(k) = alpha * xi2_m1(k) - xi2_m1(k-1);
        s3_m1(k) = alpha * xi3_m1(k) - xi3_m1(k-1);
        s4_m1(k) = alpha * xi4_m1(k) - xi4_m1(k-1);
    end

    % MFAC updates
    if k == 1
        mfa1_m1(k) = 0; mfa2_m1(k) = 0; mfa3_m1(k) = 0; mfa4_m1(k) = 0;
    else
        mfa1_m1(k) = mfa1_m1(k-1) + (rho * phi1_m1(k)) / (lamda + abs(phi1_m1(k)^2)) * xi1_m1(k);
        mfa2_m1(k) = mfa2_m1(k-1) + (rho * phi2_m1(k)) / (lamda + abs(phi2_m1(k)^2)) * xi2_m1(k);
        mfa3_m1(k) = mfa3_m1(k-1) + (rho * phi3_m1(k)) / (lamda + abs(phi3_m1(k)^2)) * xi3_m1(k);
        mfa4_m1(k) = mfa4_m1(k-1) + (rho * phi4_m1(k)) / (lamda + abs(phi4_m1(k)^2)) * xi4_m1(k);
    end

    % SMC updates
    if k == 1
        sm1_m1(k) = 0; sm2_m1(k) = 0; sm3_m1(k) = 0; sm4_m1(k) = 0;
    else
        sm1_m1(k) = sm1_m1(k-1) + (beta_m1 * phi1_m1(k)) / (sigma_m1 + (phi1_m1(k))^2) * ...
            ( (xi1_m1(k) + (y4_m1(k) - y4_m1(k-1)) + (yd(k+1) - yd(k))) / (1 + 1) ...
            - (xi1_m1(k)) / (alpha * 2) + tau_m1 * sign(s1_m1(k)) );
        sm2_m1(k) = sm2_m1(k-1) + (beta_m1 * phi2_m1(k)) / (sigma_m1 + (phi2_m1(k))^2) * ...
            ( (xi2_m1(k) + (y1_m1(k) - y1_m1(k-1)) + (y3_m1(k) - y3_m1(k-1))) / (1 + 1) ...
            - (xi2_m1(k)) / (alpha * 2) + tau_m1 * sign(s2_m1(k)) );
        sm3_m1(k) = sm3_m1(k-1) + (beta_m1 * phi3_m1(k)) / (sigma_m1 + (phi3_m1(k))^2) * ...
            ( (xi3_m1(k) + (y2_m1(k) - y2_m1(k-1)) + (yd(k+1) - yd(k))) / (1 + 1) ...
            - (xi3_m1(k)) / (alpha * 2) + tau_m1 * sign(s3_m1(k)) );
        sm4_m1(k) = sm4_m1(k-1) + (beta_m1 * phi4_m1(k)) / (sigma_m1 + (phi4_m1(k))^2) * ...
            ( (xi4_m1(k) + (y1_m1(k) - y1_m1(k-1)) + (y3_m1(k) - y3_m1(k-1))) / (1 + 1) ...
            - (xi4_m1(k)) / (alpha * 2) + tau_m1 * sign(s4_m1(k)) );
    end

    % Control signal
    if k == 1
        u1_m1(k) = 0; u2_m1(k) = 0; u3_m1(k) = 0; u4_m1(k) = 0;
    else
        u1_m1(k) = mfa1_m1(k) + gamma1_m1 * sm1_m1(k);
        u2_m1(k) = mfa2_m1(k) + gamma2_m1 * sm2_m1(k);
        u3_m1(k) = mfa3_m1(k) + gamma3_m1 * sm3_m1(k);
        u4_m1(k) = mfa4_m1(k) + gamma4_m1 * sm4_m1(k);
    end
    if k == 1
        y1_m1(k) = 0; y2_m1(k) = 0; y3_m1(k) = 0; y4_m1(k) = 0;
    end

    % Plant model update
    b1_m1 = 1.2 * n_m1 / (rT * 0.2);
    b2_m1 = 1.15 * n_m1 / (rT * 0.2);
    b3_m1 = 1.2 * n_m1 / (rT * 0.2);
    b4_m1 = 1.15 * n_m1 / (rT * 0.2);
    y1_m1(k+1) = a_m1 * y1_m1(k) + b1_m1 * u1_m1(k) - nonlinearity1 * y1_m1(k)^3 + ff_gain_m1;
    y2_m1(k+1) = a_m1 * y2_m1(k) + b2_m1 * u2_m1(k) - nonlinearity2 * y2_m1(k)^2 + ff_gain_m1;
    y3_m1(k+1) = a_m1 * y3_m1(k) + b3_m1 * u3_m1(k) - nonlinearity3 * y3_m1(k)^3 + ff_gain_m1;
    y4_m1(k+1) = a_m1 * y4_m1(k) + b4_m1 * u4_m1(k) - nonlinearity4 * y4_m1(k)^2 + ff_gain_m1;
end

% Method 2 (MFAC only) Simulation
for k = 1:m
    % Adaptive Gain update
    if k == 1
        phi1_m2(k) = 3.0; phi2_m2(k) = 3.0; phi3_m2(k) = 3.0; phi4_m2(k) = 3.0;
    elseif k == 2
        phi1_m2(k) = phi1_m2(k-1) + (eta * u1_m2(k-1) / (mu + u1_m2(k-1)^2)) * (y1_m2(k) - phi1_m2(k-1)*u1_m2(k-1));
        phi2_m2(k) = phi2_m2(k-1) + (eta * u2_m2(k-1) / (mu + u2_m2(k-1)^2)) * (y2_m2(k) - phi2_m2(k-1)*u2_m2(k-1));
        phi3_m2(k) = phi3_m2(k-1) + (eta * u3_m2(k-1) / (mu + u3_m2(k-1)^2)) * (y3_m2(k) - phi3_m2(k-1)*u3_m2(k-1));
        phi4_m2(k) = phi4_m2(k-1) + (eta * u4_m2(k-1) / (mu + u4_m2(k-1)^2)) * (y4_m2(k) - phi4_m2(k-1)*u4_m2(k-1));
    else
        phi1_m2(k) = phi1_m2(k-1) + (eta * (u1_m2(k-1) - u1_m2(k-2)) / (mu + (u1_m2(k-1) - u1_m2(k-2))^2)) * (y1_m2(k) - y1_m2(k-1) - phi1_m2(k-1) * (u1_m2(k-1) - u1_m2(k-2)));
        phi2_m2(k) = phi2_m2(k-1) + (eta * (u2_m2(k-1) - u2_m2(k-2)) / (mu + (u2_m2(k-1) - u2_m2(k-2))^2)) * (y2_m2(k) - y2_m2(k-1) - phi2_m2(k-1) * (u2_m2(k-1) - u2_m2(k-2)));
        phi3_m2(k) = phi3_m2(k-1) + (eta * (u3_m2(k-1) - u3_m2(k-2)) / (mu + (u3_m2(k-1) - u3_m2(k-2))^2)) * (y3_m2(k) - y3_m2(k-1) - phi3_m2(k-1) * (u3_m2(k-1) - u3_m2(k-2)));
        phi4_m2(k) = phi4_m2(k-1) + (eta * (u4_m2(k-1) - u4_m2(k-2)) / (mu + (u4_m2(k-1) - u4_m2(k-2))^2)) * (y4_m2(k) - y4_m2(k-1) - phi4_m2(k-1) * (u4_m2(k-1) - u4_m2(k-2)));
    end

    % Stability protection
    if k > 2 && (abs(phi1_m2(k)) <= epsilon || abs(u1_m2(k-1) - u1_m2(k-2)) <= epsilon || sign(phi1_m2(k)) ~= sign(phi1_m2(1)))
        phi1_m2(k) = phi1_m2(1);
    end
    if k > 2 && (abs(phi2_m2(k)) <= epsilon || abs(u2_m2(k-1) - u2_m2(k-2)) <= epsilon || sign(phi2_m2(k)) ~= sign(phi2_m2(1)))
        phi2_m2(k) = phi2_m2(1);
    end
    if k > 2 && (abs(phi3_m2(k)) <= epsilon || abs(u3_m2(k-1) - u3_m2(k-2)) <= epsilon || sign(phi3_m2(k)) ~= sign(phi3_m2(1)))
        phi3_m2(k) = phi3_m2(1);
    end
    if k > 2 && (abs(phi4_m2(k)) <= epsilon || abs(u4_m2(k-1) - u4_m2(k-2)) <= epsilon || sign(phi4_m2(k)) ~= sign(phi4_m2(1)))
        phi4_m2(k) = phi4_m2(1);
    end

    % Error dynamics
    xi1_m2(k) = yd(k) - 2*y1_m2(k) + y4_m2(k);
    xi2_m2(k) = y1_m2(k) - 2*y2_m2(k) + y3_m2(k);
    xi3_m2(k) = y2_m2(k) + yd(k) - 2*y3_m2(k);
    xi4_m2(k) = y1_m2(k) + y3_m2(k) - 2*y4_m2(k);

    % MFAC updates
    if k == 1
        mfa1_m2(k) = 0; mfa2_m2(k) = 0; mfa3_m2(k) = 0; mfa4_m2(k) = 0;
    else
        mfa1_m2(k) = mfa1_m2(k-1) + (rho * phi1_m2(k)) / (lamda + abs(phi1_m2(k)^2)) * xi1_m2(k);
        mfa2_m2(k) = mfa2_m2(k-1) + (rho * phi2_m2(k)) / (lamda + abs(phi2_m2(k)^2)) * xi2_m2(k);
        mfa3_m2(k) = mfa3_m2(k-1) + (rho * phi3_m2(k)) / (lamda + abs(phi3_m2(k)^2)) * xi3_m2(k);
        mfa4_m2(k) = mfa4_m2(k-1) + (rho * phi4_m2(k)) / (lamda + abs(phi4_m2(k)^2)) * xi4_m2(k);
    end

    % Control signal
    if k == 1
        u1_m2(k) = 0; u2_m2(k) = 0; u3_m2(k) = 0; u4_m2(k) = 0;
    else
        u1_m2(k) = mfa1_m2(k);
        u2_m2(k) = mfa2_m2(k);
        u3_m2(k) = mfa3_m2(k);
        u4_m2(k) = mfa4_m2(k);
    end
    if k == 1
        y1_m2(k) = 0; y2_m2(k) = 0; y3_m2(k) = 0; y4_m2(k) = 0;
    end

    % Plant model update
    b1_m2 = 1.2 * n_m2 / (rT * 0.2);
    b2_m2 = 1.15 * n_m2 / (rT * 0.2);
    b3_m2 = 1.2 * n_m2 / (rT * 0.2);
    b4_m2 = 1.15 * n_m2 / (rT * 0.2);
    y1_m2(k+1) = a_m2 * y1_m2(k) + b1_m2 * u1_m2(k) - nonlinearity1 * y1_m2(k)^3 + ff_gain_m2;
    y2_m2(k+1) = a_m2 * y2_m2(k) + b2_m2 * u2_m2(k) - nonlinearity2 * y2_m2(k)^2 + ff_gain_m2;
    y3_m2(k+1) = a_m2 * y3_m2(k) + b3_m2 * u3_m2(k) - nonlinearity3 * y3_m2(k)^3 + ff_gain_m2;
    y4_m2(k+1) = a_m2 * y4_m2(k) + b4_m2 * u4_m2(k) - nonlinearity4 * y4_m2(k)^2 + ff_gain_m2;
end

% Calculate MSE for both methods
mse_xi1_m1 = mean(xi1_m1.^2); mse_xi2_m1 = mean(xi2_m1.^2);
mse_xi3_m1 = mean(xi3_m1.^2); mse_xi4_m1 = mean(xi4_m1.^2);
mse_xi1_m2 = mean(xi1_m2.^2); mse_xi2_m2 = mean(xi2_m2.^2);
mse_xi3_m2 = mean(xi3_m2.^2); mse_xi4_m2 = mean(xi4_m2.^2);

% Print MSE values
fprintf('Method 1 (MFAC + SMC) MSE:\n');
fprintf('xi1: %.10e\n', mse_xi1_m1);
fprintf('xi2: %.10e\n', mse_xi2_m1);
fprintf('xi3: %.10e\n', mse_xi3_m1);
fprintf('xi4: %.10e\n\n', mse_xi4_m1);
fprintf('Method 2 (MFAC only) MSE:\n');
fprintf('xi1: %.10e\n', mse_xi1_m2);
fprintf('xi2: %.10e\n', mse_xi2_m2);
fprintf('xi3: %.10e\n', mse_xi3_m2);
fprintf('xi4: %.10e\n\n', mse_xi4_m2);

% Time vector
t = 1:m+1; t_err = 1:m;

% Plot outputs (y1, y2, y3, y4) for both methods
figure('Position', [100, 100, 1450, 800]);
% zoom_x_start = 90.2; % Start of zoomed x-range
% zoom_x_end = 96.5; % End of zoomed x-range

zoom_x_start = 20; % Start of zoomed x-range
zoom_x_end = 50; % End of zoomed x-range
for i = 1:4
    subplot(2, 2, i);
    plot(t, yd, '--b', 'LineWidth', 3); hold on;
    if i == 1
        plot(t, y1_m2, '-.r', 'LineWidth', 3);
        plot(t, y1_m1, '--g', 'LineWidth', 3);
        
        title('Agent 1'); grid off;
        legend('y_d(k)',  'Method [1]','Proposed scheme','orientation', 'horizontal');
        set(gca, 'FontSize', font_size);
        xlim([0 m]); % X-axis starts from 0
        % ylim([-1.5 4]); % Y-axis limits for Agent 1
        ylabel('Tracking performance', 'FontSize', font_size);
        xlabel('Time step(k)', 'FontSize', font_size);


        axes('Position', [0.20,0.65,0.13,0.10]);
        box on; hold on;
        plot(t, yd, '-.b', 'LineWidth', 3)
        plot(t, y1_m2, '-.r', 'LineWidth', 3);
        plot(t, y1_m1, '--g', 'LineWidth', 3);
        xlim([zoom_x_start zoom_x_end]);
        % yticks([-0.4,0,0.2]);
        set(gca, 'FontSize', font_size);

        
        
    elseif i == 2
        plot(t, y2_m2, '-.r', 'LineWidth', 3);
        plot(t, y2_m1, '--g', 'LineWidth', 3);


        title('Agent 2'); grid off;
        legend('y_d(k)',  'Method [1]','Proposed scheme','orientation', 'horizontal');
        set(gca, 'FontSize', font_size);
        xlim([0 m]); % X-axis starts from 0
        % ylim([-1.5 4]); % Y-axis limits for Agent 1
        ylabel('Tracking performance', 'FontSize', font_size);
        xlabel('Time step(k)', 'FontSize', font_size);

        axes('Position', [0.65,0.63,0.13,0.10]);
        box on; hold on;
        plot(t, yd, '-.b', 'LineWidth', 3)
        plot(t, y2_m2, '-.r', 'LineWidth', 3);
        plot(t, y2_m1, '--g', 'LineWidth', 3);
        xlim([zoom_x_start zoom_x_end]);
        % yticks([-0.4,0,0.2]);
        set(gca, 'FontSize', font_size);
    elseif i == 3
        plot(t, y3_m2, '-.r', 'LineWidth', 3);
        plot(t, y3_m1, '--g', 'LineWidth', 3);

        title('Agent 3'); grid off;
        legend('y_d(k)',  'Method [1]','Proposed scheme','orientation', 'horizontal');
        set(gca, 'FontSize', font_size);
        xlim([0 m]); % X-axis starts from 0
        % ylim([-1.5 4]); % Y-axis limits for Agent 1
        ylabel('Tracking performance', 'FontSize', font_size);
        xlabel('Time step(k)', 'FontSize', font_size);

        axes('Position', [0.20,0.175,0.13,0.10]);
        box on; hold on;
        plot(t, yd, '-.b', 'LineWidth', 3)
        plot(t, y3_m2, '-.r', 'LineWidth', 3);
        plot(t, y3_m1, '--g', 'LineWidth', 3);
        xlim([zoom_x_start zoom_x_end]);
        % yticks([-0.4,0,0.2]);
        set(gca, 'FontSize', font_size);
    else
        plot(t, y4_m2, '-.r', 'LineWidth', 3);
        plot(t, y4_m1, '--g', 'LineWidth', 3);

        title('Agent 4'); grid off;
        legend('y_d(k)',  'Method [1]','Proposed scheme','orientation', 'horizontal');
        set(gca, 'FontSize', font_size);
        xlim([0 m]); % X-axis starts from 0
        % ylim([-1.5 4]); % Y-axis limits for Agent 1
        ylabel('Tracking performance', 'FontSize', font_size);
        xlabel('Time step(k)', 'FontSize', font_size);

        axes('Position', [0.65,0.155,0.13,0.10]);
        box on; hold on;
        plot(t, yd, '-.b', 'LineWidth', 3)
        plot(t, y4_m2, '-.r', 'LineWidth', 3);
        plot(t, y4_m1, '--g', 'LineWidth', 3);
        xlim([zoom_x_start zoom_x_end]);
        % yticks([-0.4,0,0.2]);
        set(gca, 'FontSize', font_size);
    end
    % grid off;
    % legend('y_d(k)',  'y (MFAC only)','y (MFAC + SMC)', 'Orientation', 'horizontal');
    % set(gca, 'FontSize', font_size);

end

% Plot error dynamics (xi1, xi2, xi3, xi4) for both methods
% figure('Position', [100, 100, 1500, 750]);
% for i = 1:4
%     subplot(2, 2, i);
%     if i == 1
%         plot(t_err, xi1_m1, '-.g', 'LineWidth', 2.5); hold on;
%         plot(t_err, xi1_m2, '-b', 'LineWidth', 2.5);
%         title('Error Dynamics \xi_1(k)');
%     elseif i == 2
%         plot(t_err, xi2_m1, '-.g', 'LineWidth', 2.5); hold on;
%         plot(t_err, xi2_m2, '-b', 'LineWidth', 2.5);
%         title('Error Dynamics \xi_2(k)');
%     elseif i == 3
%         plot(t_err, xi3_m1, '-.g', 'LineWidth', 2.5); hold on;
%         plot(t_err, xi3_m2, '-b', 'LineWidth', 2.5);
%         title('Error Dynamics \xi_3(k)');
%     else
%         plot(t_err, xi4_m1, '-.g', 'LineWidth', 2.5); hold on;
%         plot(t_err, xi4_m2, '-b', 'LineWidth', 2.5);
%         title('Error Dynamics \xi_4(k)');
%     end
%     grid off;
%     legend('\xi (MFAC + SMC)', '\xi (MFAC only)', 'Orientation', 'horizontal');
%     set(gca, 'FontSize', font_size);
%     xlim([0 m]);
%     ylim([-2 5]);
% end