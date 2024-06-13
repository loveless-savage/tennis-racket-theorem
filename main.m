clear; close all;
% define tennis racket

% calculate inertia tensor
I = zeros(3);
I1 = 1; I2 = 2; I3 = 3; % test case intertial moment values to see if ODE works
% define Euler equations
dwdt = @(t,w) [(I2-I3)/I1*w(2)*w(3);(I3-I1)/I2*w(3)*w(1);(I1-I2)/I1*w(1)*w(2)];
% increase precision + detect flip events
opts = odeset('RelTol',1e-6,'Events',@flipEvent);
% solve system of ODEs
[t,w,tflip,wflip,iflip] = ode45(dwdt,[0 30],[0.01 1 0],opts);

plot(t,w);
grid on;
hold on;
plot(tflip,wflip(:,2),'*');
title('Three Axes of Angular Velocity');
xlabel('t (s)','Interpreter','latex');
ylabel('$\omega$ (rad/s)','Interpreter','latex');
legend('\(\omega_x\)','\(\omega_y\)','\(\omega_z\)','Full Period', ...
    Interpreter = 'latex');

%% Animation Station
te = 0:0.1:30; % evenly spaced time array
we = interp1(t,w(:,1),te,'spline'); % fix when we know which column of w is theta
A = max(we);

tau=te(2)-te(1);
for istep=1:length(te)
    % Position of swing relative to the pivot
    theta = we(istep);
    xswing=l1*sin(theta);
    yswing=-l1*cos(theta);
    % position of head/legs with respect to swing
    phi=A + A*cos(wp*te(istep));
    xpers=l2*sin(phi + theta);
    ypers=l2*cos(phi + theta);
    % Plot the swing and the swinger
    plot([0, xswing],[L,L+yswing],...
        [xswing+xpers,xswing-xpers],[L+yswing-ypers,L+yswing+ypers])
    % Make the x and y dimensions scale equally
    axis([-L/2 L/2 0 L])
    axis square
    % We'd like the plots frames to show at intervals of tau so the movie
    % matches the physical time scale.  However, the calculations
    % and plotting take some time, so we decrease the pause a bit.
    % Depending on the speed of your computer, you may need to adjust
    % this offset some.
    pause(tau-0.01)
end
