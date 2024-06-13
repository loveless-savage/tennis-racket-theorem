clear; close all;
% define tennis racket

% calculate inertia tensor
I = zeros(3);
I1 = 1; I2 = 2; I3 = 3; % test case intertial moment values to see if ODE works
% define Euler equations
dwdt = @(t,w) [(I2-I3)/I1*w(2)*w(3);(I3-I1)/I2*w(3)*w(1);(I1-I2)/I1*w(1)*w(2)];
% solve system of ODEs
[t,w] = ode45(dwdt,[0 10],[1 0.1 0.1],odeset('RelTol',1e-6));