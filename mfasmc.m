% Parameters initialization
d = 5;
rho = 2.5;
eta = 1.5;
lambda = 5;
mu = 1;
epsilon = 10^(-5);
alpha = 500;
T = 0.1;
gamma1 = 0.45;
gamma2 = 0.45;
gamma3 = 0.45;
gamma4 = 0.45;
omega = 0.2;
sigma = 10;
tau = 0.5;

rT = 1024;
m = 1024;
L = 200;

yd = zeros(L + 1, 1);
for k = 1:L
    yd(k) = 0.5 * sin(0.07 * pi * (k-1)) + 0.7 * cos(0.04 * pi * (k-1));
end

% Initialize arrays
phi1 = zeros(L, 1);
phi2 = zeros(L, 1);
phi3 = zeros(L, 1);
phi4 = zeros(L, 1);

mfa1 = zeros(L, 1);
mfa2 = zeros(L, 1);
mfa3 = zeros(L, 1);
mfa4 = zeros(L, 1);

sm1 = zeros(L, 1);
sm2 = zeros(L, 1);
sm3 = zeros(L, 1);
sm4 = zeros(L, 1);

u1 = zeros(L, 1);
u2 = zeros(L, 1);
u3 = zeros(L, 1);
u4 = zeros(L, 1);

y1 = zeros(L + 1, 1);
y2 = zeros(L + 1, 1);
y3 = zeros(L + 1, 1);
y4 = zeros(L + 1, 1);

xi1 = zeros(L, 1);
xi2 = zeros(L, 1);
xi3 = zeros(L, 1);
xi4 = zeros(L, 1);

s1 = zeros(L, 1);
s2 = zeros(L, 1);
s3 = zeros(L, 1);
s4 = zeros(L, 1);

% Simulation loop
for k = 1:L-1
    if k == 1
        phi1(1) = 1;
        phi2(1) = 1;
        phi3(1) = 1;
        phi4(1) = 1;
    elseif k == 2
        phi1(k) = phi1(k - 1) + eta * u1(k - 1) / (mu + u1(k - 1)^2) * (y1(k) - y1(k - 1) - phi1(k - 1) * u1(k - 1));
        phi2(k) = phi2(k - 1) + eta * u2(k - 1) / (mu + u2(k - 1)^2) * (y2(k) - y2(k - 1) - phi2(k - 1) * u2(k - 1));
        phi3(k) = phi3(k - 1) + eta * u3(k - 1) / (mu + u3(k - 1)^2) * (y3(k) - y3(k - 1) - phi3(k - 1) * u3(k - 1));
        phi4(k) = phi4(k - 1) + eta * u4(k - 1) / (mu + u4(k - 1)^2) * (y4(k) - y4(k - 1) - phi4(k - 1) * u4(k - 1));
    else
        phi1(k) = phi1(k - 1) + eta * (u1(k - 1) - u1(k - 2)) / (mu + (u1(k - 1) - u1(k - 2))^2) * ...
                  (y1(k) - y1(k - 1) - phi1(k - 1) * (u1(k - 1) - u1(k - 2)));
        phi2(k) = phi2(k - 1) + eta * (u2(k - 1) - u2(k - 2)) / (mu + (u2(k - 1) - u2(k - 2))^2) * ...
                  (y2(k) - y2(k - 1) - phi2(k - 1) * (u2(k - 1) - u2(k - 2)));
        phi3(k) = phi3(k - 1) + eta * (u3(k - 1) - u3(k - 2)) / (mu + (u3(k - 1) - u3(k - 2))^2) * ...
                  (y3(k) - y3(k - 1) - phi3(k - 1) * (u3(k - 1) - u3(k - 2)));
        phi4(k) = phi4(k - 1) + eta * (u4(k - 1) - u4(k - 2)) / (mu + (u4(k - 1) - u4(k - 2))^2) * ...
                  (y4(k) - y4(k - 1) - phi4(k - 1) * (u4(k - 1) - u4(k - 2)));
    end
    
    
    % Stability checks
    if k > 1 && (abs(phi1(k)) <= epsilon || abs(u1(k - 1) - u1(k - 2)) <= epsilon || sign(phi1(k)) ~= sign(phi1(1)))
        phi1(k) = phi1(1);
    end
    
    if k > 1 && (abs(phi2(k)) <= epsilon || abs(u2(k - 1) - u2(k - 2)) <= epsilon || sign(phi2(k)) ~= sign(phi2(1)))
        phi2(k) = phi2(1);
    end
    
    if k > 1 && (abs(phi3(k)) <= epsilon || abs(u3(k - 1) - u3(k - 2)) <= epsilon || sign(phi3(k)) ~= sign(phi3(1)))
        phi3(k) = phi3(1);
    end
    
    if k > 1 && (abs(phi4(k)) <= epsilon || abs(u4(k - 1) - u4(k - 2)) <= epsilon || sign(phi4(k)) ~= sign(phi4(1)))
        phi4(k) = phi4(1);
    end
    
    % Error calculation
    xi1(k) = yd(k) - 2 * y1(k) + y4(k);
    xi2(k) = y1(k) - 2 * y2(k) + y3(k);
    xi3(k) = y2(k) + yd(k) - 2 * y3(k);
    xi4(k) = y1(k) + y3(k) - 2 * y4(k);
    
    % Stability term
    s1(k) = alpha * xi1(k) - xi1(k-1);
    s2(k) = alpha * xi2(k) - xi2(k-1);
    s3(k) = alpha * xi3(k) - xi3(k-1);
    s4(k) = alpha * xi4(k) - xi4(k-1);
    
    % MFA calculation
    mfa1(k) = mfa1(k - 1) + (rho * phi1(k)) / (lambda + abs(phi1(k))^2) * xi1(k);
    mfa2(k) = mfa2(k - 1) + (rho * phi2(k)) / (lambda + abs(phi2(k))^2) * xi2(k);
    mfa3(k) = mfa3(k - 1) + (rho * phi3(k)) / (lambda + abs(phi3(k))^2) * xi3(k);
    mfa4(k) = mfa4(k - 1) + (rho * phi4(k)) / (lambda + abs(phi4(k))^2) * xi4(k);
    
    % SM calculation
    sm1(k) = sm1(k - 1) + (omega * phi1(k)) / (sigma + abs(phi1(k))^2) * ((alpha * (xi1(k)) - xi1(k)) / (alpha * 2) - y1(k) + tau * sign(s1(k)));
    sm2(k) = sm2(k - 1) + (omega * phi2(k)) / (sigma + abs(phi2(k))^2) * ((alpha * (xi2(k)) - xi2(k)) / (alpha * 2) - y2(k) + tau * sign(s2(k)));
    sm3(k) = sm3(k - 1) + (omega * phi3(k)) / (sigma + abs(phi3(k))^2) * ((alpha * (xi3(k)) - xi3(k)) / (alpha * 2) - y3(k) + tau * sign(s3(k)));
    sm4(k) = sm4(k - 1) + (omega * phi4(k)) / (sigma + abs(phi4(k))^2) * ((alpha * (xi4(k)) - xi4(k)) / (alpha * 2) - y4(k) + tau * sign(s4(k)));
    
    % Update u values
    u1(k) = mfa1(k) + gamma1 * sm1(k);
    u2(k) = mfa2(k) + gamma2 * sm2(k);
    u3(k) = mfa3(k) + gamma3 * sm3(k);
    u4(k) = mfa4(k) + gamma4 * sm4(k);
    
    % Update outputs
    y1(k + 1) = m / (rT * 1) * u1(k);
    y2(k + 1) = m / (rT * 1) * u2(k);
    y3(k + 1) = m / (rT * 1) * u3(k);
    y4(k + 1) = m / (rT * 1) * u4(k);
end

% Plotting
figure; 
hold on;
plot(yd(1:end-1), '-y', 'LineWidth', 2, 'DisplayName', '$y_d$');
plot(y1(1:end-1), '--r', 'LineWidth', 2, 'DisplayName', '$y_1$');
plot(y2(1:end-1), '-.b', 'LineWidth', 2, 'DisplayName', '$y_2$');
plot(y3(1:end-1), '--k', 'LineWidth', 2, 'DisplayName', '$y_3$');
plot(y4(1:end-1), '-.g', 'LineWidth', 2, 'DisplayName', '$y_4$');
hold off;

% Labels and legend
xlabel('Time step', 'FontSize', 14);
ylabel('Output', 'FontSize', 14);
legend('Location', 'Best', 'FontSize', 14);
set(gca, 'FontSize', 12);
xlim([0 L]);
grid off;
tight_layout;
