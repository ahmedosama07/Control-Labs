speed_ctrl = tf([100 40], [1 0]);               % Speed Controller
trq_ctrl = tf([10 6], [1 0]);                   % Torque Controller
motor_inertia = tf(1, [7.226 0]);               % Motor Inertia
K_cs = 0.5;                                     % Current Sensor Senstivity
K_ss = 0.0433;                                  % Speed Sensor Senstivity
k_f = 0.1;                                      % Viscous friction
k_b = 2;                                        % Back EMF constant
R_a = 1;                                        % Motor Resistance
motor_trq = 1.8;                                % Motor torque
vehicle_trq = 0.6154;                           % Vehicle torque
vehicle_speed = 0.0615;                         % Vehile speed

stage1 = feedback(motor_inertia, k_f);
k_b = k_b / vehicle_speed;
K_ss = K_ss / vehicle_speed;
stage2 = series(stage1, vehicle_speed);
stage3 = feedback(stage2, vehicle_trq);
stage4 = series(motor_trq, stage3);
K_cs = K_cs / stage4;
stage5 = series((1/R_a), stage4);
stage6 = feedback(stage5, k_b);
stage7 = series(trq_ctrl, stage6);
stage8 = feedback(stage7, K_cs);
stage9 = series(speed_ctrl, stage8);
hev = feedback(stage9, K_ss);

[num, den] = tfdata(hev);
transfere_function = tf(num{1}/den{1}(1), den{1}/den{1}(1));
display(transfere_function);