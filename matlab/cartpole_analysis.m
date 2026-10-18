%% Cart-pole: nonlinear model, linearization, stability and controllability
% A pendulum on a motorized cart, with state X = [x; v; theta; omega].
% theta is measured from the downward (hanging) position: theta = 0 hangs
% straight down and theta = pi is the inverted (upright) position.
% The only input is a horizontal force u on the cart.

clear; clc; close all;

%% Parameters
M = 4;          % cart mass (kg)
m = 0.6;        % pendulum mass (kg)
L = 1.2;        % pendulum arm length (m)
d = 0.4;        % viscous friction (damping) coefficient of the cart
g = 9.81;       % gravitational acceleration (m/s^2)
u = 0;          % free response (no control force)

% Sign convention: in the pendcart dynamics below, g enters with the sign
% used by the standard cart-pole demo this model follows (g negative), so
% the calls pass -g.

%% Nonlinear model: free-response simulations
% Case 1: released from rest 0.5 rad from the downward position
% Case 2: released from rest from the horizontal position
% Case 3: released from rest 0.2 rad from the upright position
% Case 4: released at 0.5 rad with an initial angular velocity of 3 rad/s
X0cases = [0, 0, 0.5,     0;
           0, 0, pi/2,    0;
           0, 0, pi - 0.2, 0;
           0, 0, 0.5,     3];
Tcases = [12, 12, 10, 12];

for k = 1:4
    X0 = X0cases(k,:);
    tspan = 0:0.01:Tcases(k);
    [t, X] = ode45(@(t,X) pendcart(X, m, M, L, -g, d, u), tspan, X0);

    figure;
    subplot(1,2,1);
    plot(t, rad2deg(X(:,3)), 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('\theta (deg)');
    title('Pendulum angle \theta(t) (0 = straight down)');
    subplot(1,2,2);
    plot(t, X(:,1), 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('x (m)');
    title('Cart position x(t)');
    sgtitle(sprintf('Nonlinear model, Case %d', k));
end

%% Linearized models, eigenvalues and stability
% Linearizing about the two equilibria (deviation variables are used for
% the upper equilibrium, theta = pi):

A_lower = [0,       1,      0,             0;
           0, -d/M,    m*g/M,           0;
           0,       0,      0,           1;
           0,  d/(M*L), -(m+M)*g/(M*L),  0];

A_upper = [0,       1,      0,             0;
           0, -d/M,    m*g/M,           0;
           0,       0,      0,           1;
           0, -d/(M*L),  (m+M)*g/(M*L),  0];

disp('Eigenvalues, lower equilibrium:');
disp(eig(A_lower));
disp('Eigenvalues, upper equilibrium:');
disp(eig(A_upper));

%% Linearized vs. nonlinear response

% Lower equilibrium: initial angle deviation of 0.5 rad (as in Case 1)
tspan = 0:0.01:12;
[tNL, XNL] = ode45(@(t,X) pendcart(X, m, M, L, -g, d, u), tspan, [0;0;0.5;0]);
[tLI, ZLI] = ode45(@(t,Z) A_lower*Z, tspan, [0;0;0.5;0]);

figure;
subplot(1,2,1);
plot(tNL, rad2deg(XNL(:,3)), 'LineWidth', 1.5); hold on;
plot(tLI, rad2deg(ZLI(:,3)), '--', 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('\theta (deg)');
title('Pendulum angle (0 = straight down)');
legend('Nonlinear', 'Linearized (lower)', 'Location', 'best');
subplot(1,2,2);
plot(tNL, XNL(:,1), 'LineWidth', 1.5); hold on;
plot(tLI, ZLI(:,1), '--', 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('x (m)');
title('Cart position');
legend('Nonlinear', 'Linearized (lower)', 'Location', 'best');
sgtitle('Lower equilibrium: nonlinear vs. linearized model');

% Upper equilibrium: initial angle deviation of -0.2 rad from upright
tspan = 0:0.01:1;
[tNU, XNU] = ode45(@(t,X) pendcart(X, m, M, L, -g, d, u), tspan, [0;0;pi-0.2;0]);
[tLU, ZLU] = ode45(@(t,Z) A_upper*Z, tspan, [0;0;-0.2;0]);

figure;
subplot(1,2,1);
plot(tNU, rad2deg(XNU(:,3) - pi), 'LineWidth', 1.5); hold on;
plot(tLU, rad2deg(ZLU(:,3)), '--', 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('\theta - \pi (deg)');
title('Angle deviation from the upright');
legend('Nonlinear', 'Linearized (upper)', 'Location', 'best');
subplot(1,2,2);
plot(tNU, XNU(:,1), 'LineWidth', 1.5); hold on;
plot(tLU, ZLU(:,1), '--', 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('x (m)');
title('Cart position');
legend('Nonlinear', 'Linearized (upper)', 'Location', 'best');
sgtitle('Upper equilibrium: nonlinear vs. linearized model');

% Phase plane about the lower equilibrium, theta(0) = 0.9 rad
tspan = 0:0.01:10;
[tNP, XNP] = ode45(@(t,X) pendcart(X, m, M, L, -g, d, u), tspan, [0;0;0.9;0]);
[tLP, ZLP] = ode45(@(t,Z) A_lower*Z, tspan, [0;0;0.9;0]);

figure;
plot(rad2deg(XNP(:,3)), rad2deg(XNP(:,4)), 'LineWidth', 1.5); hold on;
plot(rad2deg(ZLP(:,3)), rad2deg(ZLP(:,4)), '--', 'LineWidth', 1.5);
grid on; xlabel('Pendulum angle \theta (deg)');
ylabel('Angular velocity \omega (deg/s)');
title('Phase plane about the lower equilibrium');
legend('Nonlinear', 'Linearized (lower)', 'Location', 'best');

%% Controllability
% The input enters through the cart force only. Linearizing the input
% term of the nonlinear model about each equilibrium gives:

B_lower = [0; 1/M; 0; -1/(M*L)];
B_upper = [0; 1/M; 0;  1/(M*L)];

C_lower = ctrb(A_lower, B_lower);   % C = [B, A*B, A^2*B, A^3*B]
C_upper = ctrb(A_upper, B_upper);
fprintf('Lower equilibrium:  rank(C) = %d, det(C) = %.4f\n', ...
        rank(C_lower), det(C_lower));
fprintf('Upper equilibrium:  rank(C) = %d, det(C) = %.4f\n', ...
        rank(C_upper), det(C_upper));
% Both controllability matrices have rank 4, so both linearized
% systems are controllable from the single cart-force input.

%% Nonlinear dynamics (pendcart)
function dx = pendcart(x,m,M,L,g,d,u)

Sx = sin(x(3));
Cx = cos(x(3));
D = m*L*L*(M+m*(1-Cx^2));

dx(1,1) = x(2);
dx(2,1) = (1/D)*(-m^2*L^2*g*Cx*Sx + m*L^2*(m*L*x(4)^2*Sx - d*x(2))) + m*L*L*(1/D)*u;
dx(3,1) = x(4);
dx(4,1) = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*x(4)^2*Sx - d*x(2))) - m*L*Cx*(1/D)*u;
end
