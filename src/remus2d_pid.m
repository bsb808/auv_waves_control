function del = remus2d_pid(z_setpoint, state)

% Simple classical controller (PID) to illustrate the interface
% 
% Inputs:
% * z_setpoint = goal position in z direction [m]
% * state = state vector [x, z, theta, u, w, q]
% 
% Outputs:
% * del = fin deflection [rad]
% 

% Gains
Kp = 0.15; 
Kd = 0;
error = z_setpoint-state(2);
del = -1*(Kp*error + Kd*(-state(5)));

return;