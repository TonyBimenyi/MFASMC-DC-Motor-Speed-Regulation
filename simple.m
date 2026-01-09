clear; clc;

% Parameters
d = 5;
rho = 1;
eta = 1;
lambda = 20;
mu = 1;
epsilon = 1e-5;
alpha = 1;
T = 0.1;

gamma1 = 0.15;
gamma2 = 0.15;
gamma3 = 0.45;
gamma4 = 0.45;

rT = 1024;
m = 350;
L = 200;

% Generate a more complex desired trajectory
yd = zeros(L+1, 1);
for k = 1:L
    yd(k) = 0.5 * sin(k * pi / 30) + 0.3 * cos(k * pi / 10);
end

% Initialize arrays
phi1 = zeros(L,1); phi2 = zeros(L,1);
phi3 = zeros(L,1); phi4 = zeros(L,1);

mfa1 = zeros(L,1); mfa2 = zeros(L,1);
mfa3 = zeros(L,1); mfa4 = zeros(L,1);

sm1 = zeros(L,1); sm2 = zeros(L,1);
sm3 = zeros(L,1); sm4 = zeros(L,1);

u1 = zeros(L,1); u2 = zeros(L,1);
u3 = zeros(L,1); u4 = zeros(L,1);

y1 = zeros(L+1,1); y2 = zeros(L+1,1);
y3 = zeros(L+1,1); y4 = zeros(L+1,1);

e1 = zeros(L+1,1); e2 = zeros(L+1,1);
e3 = zeros(L+1,1); e4 = zeros(L+1,1);

si1 = zeros(L,1); si2 = zeros(L,1);
si3 = zeros(L,1); si4 = zeros(L,1);

% Initial outputs
y1(1) = 0.1;
y2(1) = 0.1;
y3(1) = 0.1;
y4(1) = 0.1;

% Simulation Loop
for k = 2 : L-1

    if k == 2
        phi1(k) = 1;
        phi2(k) = 1;
        phi3(k) = 1;
        phi4(k) = 1;

    elseif k == 3
        phi1(k) = phi1(k-1) + eta * u1(k-1)/(mu + u1(k-1)^2) * (y1(k)-y1(k-1)-phi1(k-1)*u1(k-1));
        phi2(k) = phi2(k-1) + eta * u2(k-1)/(mu + u2(k-1)^2) * (y2(k)-y2(k-1)-phi2(k-1)*u2(k-1));
        phi3(k) = phi3(k-1) + eta * u3(k-1)/(mu + u3(k-1)^2) * (y3(k)-y3(k-1)-phi3(k-1)*u3(k-1));
        phi4(k) = phi4(k-1) + eta * u4(k-1)/(mu + u4(k-1)^2) * (y4(k)-y4(k-1)-phi4(k-1)*u4(k-1));

    else
        phi1(k) = phi1(k-1) + (eta*(u1(k-1)-u1(k-2))/(mu + abs(u1(k-1)-u1(k-2))^2)) * ...
            (y1(k)-y1(k-1)-phi1(k-1)*(u1(k-1)-u1(k-2)));

        phi2(k) = phi2(k-1) + (eta*(u2(k-1)-u2(k-2))/(mu + abs(u2(k-1)-u2(k-2))^2)) * ...
            (y2(k)-y2(k-1)-phi2(k-1)*(u2(k-1)-u2(k-2)));

        phi3(k) = phi3(k-1) + (eta*(u3(k-1)-u3(k-2))/(mu + abs(u3(k-1)-u3(k-2))^2)) * ...
            (y3(k)-y3(k-1)-phi3(k-1)*(u3(k-1)-u3(k-2)));

        phi4(k) = phi4(k-1) + (eta*(u4(k-1)-u4(k-2))/(mu + abs(u4(k-1)-u4(k-2))^2)) * ...
            (y4(k)-y4(k-1)-phi4(k-1)*(u4(k-1)-u4(k-2)));
    end

    % Sliding errors
    si1(k) = yd(k) - 2*y1(k) + y4(k);
    si2(k) = y1(k) - 2*y2(k) + y3(k);
    si3(k) = y2(k) + yd(k) - 2*y3(k);
    si4(k) = y1(k) + y3(k) - 2*y4(k);

    % MFAC learning law
    if k > 2
        mfa1(k) = mfa1(k-1) + (rho*phi1(k))/(lambda + abs(phi1(k))^2) * si1(k);
        mfa2(k) = mfa2(k-1) + (rho*phi2(k))/(lambda + abs(phi2(k))^2) * si2(k);
        mfa3(k) = mfa3(k-1) + (rho*phi3(k))/(lambda + abs(phi3(k))^2) * si3(k);
        mfa4(k) = mfa4(k-1) + (rho*phi4(k))/(lambda + abs(phi4(k))^2) * si4(k);
    end

    % Slidinzg mode update
    if k > 2
        sm1(k) = sm1(k-1) + (yd(k+1)-y1(k) + alpha*yd(k+1) - y1(k) + epsilon*T*sign(k));
        sm2(k) = sm2(k-1) + (yd(k+1)-y1(k) + alpha*yd(k+1) - y1(k) + epsilon*T*sign(k));
        sm3(k) = sm3(k-1) + (yd(k+1)-y1(k) + alpha*yd(k+1) - y1(k) + epsilon*T*sign(k));
        sm4(k) = sm4(k-1) + (yd(k+1)-y1(k) + alpha*yd(k+1) - y1(k) + epsilon*T*sign(k));
    end

    % Control input
    if k == 2
        u1(k) = 0.1; u2(k) = 0.1; u3(k) = 0.1; u4(k) = 0.1;
    else
        u1(k) = mfa1(k) + gamma1 * sm1(k);
        u2(k) = mfa2(k) + gamma2 * sm2(k);
        u3(k) = mfa3(k) + gamma3 * sm3(k);
        u4(k) = mfa4(k) + gamma4 * sm4(k);
    end

    % Output update
    y1(k+1) = m / (rT*0.1) * u1(k);
    y2(k+1) = m / (rT*0.1) * u2(k);
    y3(k+1) = m / (rT*0.3) * u3(k);
    y4(k+1) = m / (rT*0.3) * u4(k);

end

% ========= PLOT ==========
figure;
hold on;
plot(yd, '-r', 'LineWidth',2);
plot(y1, '-r', 'LineWidth',2);
plot(y2, '-.g', 'LineWidth',2);
plot(y3, '-y', 'LineWidth',2);
plot(y4, '--k', 'LineWidth',2);

xlabel('Time Step','FontSize',14);
ylabel('Output','FontSize',14);

legend('\phi_1','\phi_2','\phi_3','\phi_4','FontSize',14);
grid on;
xlim([1 L]);
