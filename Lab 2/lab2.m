% MATLAB Script for Control System Analysis
% Author: Ahmed Osama
% Date: 26 / 11 / 2024
% Description:
% This script analyzes a control system with a second-order transfer
% function. It performs tasks such as deriving the state-space
% representation, evaluating stability, and visualizing responses for
% various system parameters.

%% Clear Workspace and Initialize
clc; clear all;  % Clear command window and variables
%% Evaluate the closed loop transfer function
% Define system parameters for a second-order control system.
J = 600000;  % Moment of inertia (kg*m^2)
B = 20000;   % Damping coefficient (N*m*s)
syms K s     % Symbolic variables for gain (K) and Laplace variable (s)

% Define the closed-loop transfer function H(s)
H(s) = K / (J * s^2 + B * s + K);

%% State-Space Representation for K = 1
% Substitute K = 1 into the transfer function and find state-space representation.
[nm, dn] = numden(subs(H, K, 1));  % Get numerator and denominator
nmd = sym2poly(nm);  % Convert symbolic numerator to polynomial coefficients
dnd = sym2poly(dn);  % Convert symbolic denominator to polynomial coefficients

tf_K1 = tf(nmd, dnd);  % Transfer function for K = 1
ss_K1 = ss(tf_K1);     % State-space representation for K = 1


% Plot root locus for K = 1
figure;
rlocusplot(tf_K1);
title('Root Locus for K = 1');

%% Maximum Value of K for Stability
% Calculate the maximum value of K that keeps the system stable.
syms Kmax s_Kmax
Kmax = -J * s_Kmax^2 - B * s_Kmax;  % Stability criterion
s_Kmax = solve(diff(Kmax, s_Kmax) == 0);  % Solve for critical s
sln = double(subs(Kmax));  % Substitute and evaluate Kmax

% Generate transfer function for max K
[nm, dn] = numden(subs(H, K, sln));
nmd = sym2poly(nm);
dnd = sym2poly(dn);
tf_Kmax = tf(nmd, dnd);

% Plot root locus for max K
figure;
rlocusplot(tf_Kmax);
title('Root Locus for Maximum Stable K');

%% Maximum K for Overshoot < 10%
% Determine K such that the system overshoot M_p is less than 10%.
syms K_Mp10 real
ineq = solve(((B^2)/(4*J) < K_Mp10) & ...
    (K_Mp10 < ((2/B) * sqrt((J * log(0.1)^2) / ((pi^2) + (log(0.1)^2)))) ^ -2), ...
    K_Mp10, 'ReturnConditions', true);

% Plot feasible K values on a number line
x = 0 : 0.01 : 1000;
figure;
hold on;
range = x(boolean(subs(ineq.conditions, x)));

% Plot base line and highlight feasible K values
plot(x, zeros(size(x)), 'k', 'LineWidth', 1);
plot(range, zeros(size(range)), 'r', 'LineWidth', 3);
plot(range(1), 0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
plot(range(end), 0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
ylim([-0.5, 0.5]);
title('Values of K for {M}_{p} < 10%');
grid on;
hold off;

%% Maximum K for Rise Time < 80 Seconds
% Binary search to find maximum K such that rise time is under 80 seconds.
RiseTime_max = 80;  % Target rise time (seconds)
low = 0; high = 1000; precision = 1e-5;  % Search range and precision
k = low;
while (high - low) > precision
    mid = (low + high) / 2;
    [nm,dn] = numden(subs(H, K, mid));
     nmd = sym2poly(nm);
     dnd = sym2poly(dn);
     sys_tf = tf(nmd, dnd);
     sys_ss = ss(sys_tf);
     info = stepinfo(sys_ss);
     
     if info.RiseTime < RiseTime_max
         k = mid;
         high = mid;
     else 
         low = mid;
     end
end

% Plot feasible K values on a number line
syms K_tr80 real
ineq = solve(((B^2)/(4*J) < K_tr80) &(K_tr80 > k), K_tr80, 'ReturnConditions',true);
x = 0 : precision : 1000;
%range = x(boolean(subs(ineq.conditions, x)));
range = x(boolean(x>k));
figure;
hold on;
plot(x, zeros(size(x)), 'k', 'LineWidth', 1);
plot(range, zeros(size(range)), 'r', 'LineWidth', 3);
plot(range(1), 0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
ylim([-0.5, 0.5]);
title('Values of K for Rise Time < 80 sec');
grid on;
hold off;
%% Step Response for Various K Values
% Analyze step response for K = 200, 400, 1000, 2000.
K_values = [200, 400, 1000, 2000];
state_spaces = [];
figure;
for i = 1:length(K_values)
    [nm,dn] = numden(subs(H, K, K_values(i)));
    nmd = sym2poly(nm);
    dnd = sym2poly(dn);
    curr_tf = tf(nmd, dnd);
    state_spaces = [state_spaces ss(curr_tf)];
    
    ss_err = 1 / (1+dcgain(curr_tf));  % Steady-state error
    info = stepinfo(state_spaces(i));  % Extract system info
    disp(['For K = ', num2str(K_values(i)), ':']);
    disp([' Overshoot: ', num2str(info.Overshoot)]);
    disp([' Rise Time: ', num2str(info.RiseTime)]);
    disp([' Steady-State Error: ', num2str(ss_err)]);
    disp(' ');
    subplot(2, 2, i);
    step(state_spaces(i));
    title(['Step Response for K = ', num2str(K_values(i))]);
end
sgtitle('Step Response for Different Values of K');

%% Pole-Zero Maps for K Values
% Visualize pole-zero locations for various K values.
figure;
for i = 1:length(K_values)
    subplot(2, 2, i);
    pzmap(state_spaces(i));
    title(['Pole-Zero Map for K = ', num2str(K_values(i))]);
end
sgtitle('Pole-Zero Map for Different Values of K');

%% Using Simulink, adding poles at -200, -20, -10, and -2
poles_step_response = load('poles_step_response.mat');
sys_info = stepinfo(poles_step_response.LinearAnalysisToolProject.Results.Data.Value);
disp('For higher-order pole is nonexistent');
disp([' Overshoot: ', num2str(sys_info(1).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(1).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(1).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(1).RiseTime)]);

disp('For higher-order pole at -200');
disp([' Overshoot: ', num2str(sys_info(2).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(2).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(2).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(2).RiseTime)]);

disp('For higher-order pole at -20');
disp([' Overshoot: ', num2str(sys_info(3).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(3).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(3).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(3).RiseTime)]);

disp('For higher-order pole at -10');
disp([' Overshoot: ', num2str(sys_info(4).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(4).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(4).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(4).RiseTime)]);

disp('For higher-order pole at -2');
disp([' Overshoot: ', num2str(sys_info(5).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(5).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(5).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(5).RiseTime)]);

%% Using Simulink, adding zeros at -200, -50, -20, -10, -5, and -2
zeros_step_response = load('zeros_step_response.mat');
sys_info = stepinfo(zeros_step_response.LinearAnalysisToolProject.Results.Data.Value);
disp('For zero is nonexistent');
disp([' Overshoot: ', num2str(sys_info(1).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(1).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(1).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(1).RiseTime)]);

disp('For zero at -200');
disp([' Overshoot: ', num2str(sys_info(2).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(2).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(2).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(2).RiseTime)]);

disp('For zero at -50');
disp([' Overshoot: ', num2str(sys_info(3).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(3).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(3).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(3).RiseTime)]);

disp('For zero at -20');
disp([' Overshoot: ', num2str(sys_info(4).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(4).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(4).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(4).RiseTime)]);

disp('For zero at -10');
disp([' Overshoot: ', num2str(sys_info(5).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(5).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(5).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(5).RiseTime)]);

disp('For zero at -5');
disp([' Overshoot: ', num2str(sys_info(6).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(6).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(6).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(6).RiseTime)]);

disp('For zero at -2');
disp([' Overshoot: ', num2str(sys_info(7).Overshoot)]);
disp([' Settling time: ', num2str(sys_info(7).SettlingTime)]);
disp([' Peak time: ', num2str(sys_info(7).PeakTime)]);
disp([' Rise time: ', num2str(sys_info(7).RiseTime)]);