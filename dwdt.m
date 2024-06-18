%% ODE function
% parameters:
% t => scalar timestamp (not used directly)
% w => differential operand
% ├─ w(1:3) => xyz components of angular velocity
% ├─ w(4:6) => rotated x-unit vector to track total rotation
% └─ w(7:9) => rotated y-unit vector
% I1,I2,I3 => moments of inertia around each principal axis
% output:
% dw => time derivative of each component of w
function dw = dwdt(t,w,I1,I2,I3)

% initialize derivative value
dw = 0*w;

% Euler's rotation equations
dw(1) = (I2-I3)/I1*w(2)*w(3);
dw(2) = (I3-I1)/I2*w(3)*w(1);
dw(3) = (I1-I2)/I1*w(1)*w(2);

% integrate basis vector motion to track total rotation
% df/dt = w x f, where f is a basis vector in a rotation matrix
dw(4:6) = cross(w(1:3),w(4:6)); % x-axis
dw(7:9) = cross(w(1:3),w(7:9)); % y-axis
% z-axis can be calculated from the first two