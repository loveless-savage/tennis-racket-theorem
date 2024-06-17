clear; close all;
% define inertia tensor
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

%% Animation
% this is definitely not perfect yet
close all;
% draw ellipse
theta = linspace(0, 2*pi, 100); % array to build ellipse over
z = I1 * ones(size(theta)); % Constant initial z-coordinate for the ellipse
x = I2 * cos(theta);
y = I3 * sin(theta);
% draw ellipse
s = fill3(x,y,z,'b');
boxlim = 4;
light('Position', [4 4 4],'Style','Infinite')
% animate ellipse spinning
for i = 1:length(w)
    % rotate can only take a 1x3 array, so extract diagonals - maybe change
    % which values we use to look at the rotation?
    spin = [rot(1,1,i),rot(2,2,i),rot(3,3,i)];
    rotate(s,spin,1);
    drawnow;
    pause(0.001);
    axis([-boxlim boxlim -boxlim boxlim -boxlim boxlim]);
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
