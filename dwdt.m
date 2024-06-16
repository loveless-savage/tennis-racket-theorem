%% ODE function
% parameters:
function dw = dwdt(t,w,I1,I2,I3)
% initialize derivative value
dw = 0*w;

% Euler's rotation equations
dw(1) = (I2-I3)/I1*w(2)*w(3);
dw(2) = (I3-I1)/I2*w(3)*w(1);
dw(3) = (I1-I2)/I1*w(1)*w(2);

% also integrate rotation vectors to track rotation of basis vectors:
% de/dt = w x e, where e is a basis vector in a rotation matrix
dw(4:6) = cross(w(1:3),w(4:6)); % 1st axis
dw(7:9) = cross(w(1:3),w(7:9)); % 2nd axis
% 3rd axis can be calculated from the first two