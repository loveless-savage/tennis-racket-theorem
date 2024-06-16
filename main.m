clear; close all;
% define tennis racket

% calculate inertia tensor
I = zeros(3);
I1 = 1; I2 = 2; I3 = 3; % test case intertial moment values to see if ODE works
% duration
tmax = 30;
% initial angular velocities: give most to w(2), but a tiny bit to w(1)
% Some cool initial conditions to try: 
% [0.01 1 0]
% [0.01 1 1]
% [0.01 10 1]
% [0.01 1 .1]
w_init = [0.01 1 0];

% increase precision + detect flip events for the ODE solver
opts = odeset('RelTol',1e-6,'Events',@flipEvent);

% solve system of ODEs
[t_raw,u,tflip,uflip,iflip] = ode45( ...
    @(t,u) dwdt(t,u,I1,I2,I3), ...
    [0 tmax], ...
    [w_init, 1 0 0, 0 1 0], ...
    opts);
% collect angular velocities relative to each principal axis
w_raw = u(:,1:3);
wflip = uflip(:,1:3);
% collect x & y rotation vectors so we can build a rotation matrix
rotx_raw = u(:,4:6);
roty_raw = u(:,7:9);

% interpolate raw values w/ evenly spaced timesteps
N = 2^ceil(log2(length(t_raw))); % next power of 2 above sample #
t = (0:N-1)*tmax/N;
w = interp1(t_raw,w_raw,t,'linear'); % FIXME: dim

% to build interpolated rotation matrices, interpolate both axes...
rotx = interp1(t_raw,rotx_raw,t,'linear'); % x-axis
roty = interp1(t_raw,roty_raw,t,'linear'); % y-axis
% ... and then extrapolate the z-axis
rotz = cross(rotx,roty,2);
% now combine all three axes into one structure
rot = cat(3,rotx,roty,rotz);
rot = permute(rot,[2,3,1]);
% to get the 3x3 rotation matrix of frame n, type rot(:,:,n)

%% show angular velocity components over time
plot(t,w);
grid on;
hold on;
plot(tflip,wflip(:,2),'*');
title('Three Axes of Angular Velocity');
xlabel('t (s)','Interpreter','latex');
ylabel('$\omega$ (rad/s)','Interpreter','latex');
legend('\(\omega_x\)','\(\omega_y\)','\(\omega_z\)','Full Period', ...
    Interpreter = 'latex',Location = 'best');

%% Animation Station
te = 0:0.1:tmax; % evenly spaced time array
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

%% Frequency Analysis via Fourier Transform
close all;

% Calculate the Fast Fourier Transform (FFT)
N = length(w); % Number of samples
F=fft(w); % FFT of the signal
P = abs(F).^2; % Power spectrum
tau=t(2)-t(1); % Time step
Fs=1/tau; % Sampling frequency
dv = Fs / N; % Frequency resolution
v = (0:N-1) * dv; % Frequency vector

semilogy(v, P)
xlabel('v')
ylabel('P(v)')
title('Power Spectrum')
legend('Axis1', 'Axis2', 'Axis3')
