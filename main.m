% Andrew Jones, Eben Lonsdale, Matthew Rundquist
clear; close all;

% define the principal moments of inertia
I1 = 1; I2 = 2; I3 = 3;
% total duration of the simulation
tmax = 30;
% initial angular velocities: give most to w(2), but a tiny bit to w(1)
% Some cool initial conditions to try: 
% w_init = [0.01 1 0];
% w_init = [0.01 1 1];
% w_init = [0.01 10 1];
% w_init = [0.01 1 .1];
 w_init = [0.01 1 0];

% increase precision for the ODE solver + detect flip events
opts = odeset('RelTol',1e-6,'Events',@flipEvent);

% solve system of ODEs
[t_raw,u,tflip,uflip,iflip] = ode45( ...
    @(t,u) dwdt(t,u,I1,I2,I3), ... % ODE
    [0 tmax], ... % time range
    [w_init, 1 0 0, 0 1 0], ... % initial angular velocities + unit axes
    opts); % precision + events
% collect angular velocities relative to each principal axis
w_raw = u(:,1:3);
wflip = uflip(:,1:3);
% determine the period, the time it takes for one flip
period=tflip;
% collect x & y rotation vectors so we can build a rotation matrix
rotx_raw = u(:,4:6);
roty_raw = u(:,7:9);

% interpolate raw values w/ evenly spaced timesteps
N = 2^ceil(log2(length(t_raw))); % next power of 2 above ode45's sample#
t = (0:N-1)*tmax/N;
w = interp1(t_raw,w_raw,t,'linear');

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
plot(tflip,wflip(:,2),'*'); % star the events marking each period
title('Three Axes of Angular Velocity');
xlabel('t (s)','Interpreter','latex');
ylabel('$\omega$ (rad/s)','Interpreter','latex');
legend('\(\omega_x\)','\(\omega_y\)','\(\omega_z\)','Full Period', ...
    Interpreter = 'latex',Location = 'best');

%% Animation
close all;
% draw ellipse
theta = linspace(0, 2*pi, 100);N = size(theta); % array to build ellipse over
x1 = I2 * cos(theta);
y1 = I3 * sin(theta);
z1 = I1 * zeros(N);
s1 = [x1;y1;z1];
x2 = I2*zeros(N);
y2 = I3*sin(theta);
z2 = I2*cos(theta);
s2 = [x2;y2;z2];
x3 = I2*cos(theta);
y3 = I3*zeros(N);
z3 = I2*sin(theta);
s3 = [x3;y3;z3];
hold on;
% draw ellipse
fill3(x1,y1,z1,'b');
fill3(x2,y2,z2,'r');
fill3(x3,y3,z3,'y')
title('Animation of Intermediate Axis Theorem');
xlabel('x');
ylabel('y');
zlabel('z');
legend('I2','I3','I1');
boxlim = 4;
axis([-boxlim boxlim -boxlim boxlim -boxlim boxlim]);
view([20 20 20])
hold off;
% animate ellipse spinning
for n = 1:length(w)
    clf;
    sn1 = rot(:,:,n)*s1;
    sn2 = rot(:,:,n)*s2;
    sn3 = rot(:,:,n)*s3;
    hold on;
    fill3(sn1(1,:),sn1(2,:),sn1(3,:),'b');
    fill3(sn2(1,:),sn2(2,:),sn2(3,:),'r');
    fill3(sn3(1,:),sn3(2,:),sn3(3,:),'y');
    title('Animation of Intermediate Axis Theorem');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    legend('I2','I3','I1');
    axis([-boxlim boxlim -boxlim boxlim -boxlim boxlim]);
    view([20 20 20])
    drawnow;
    hold off;
    pause(0.0001);
end
%% Frequency Analysis via Fourier Transform
close all;

% Calculate the Fast Fourier Transform (FFT) for each axis
N = length(t); % Number of samples
tau = t(2) - t(1); % Time step
Fs = 1 / tau; % Sampling frequency

% Initialize power spectrum matrix
P = zeros(N, 3);
v = linspace(-Fs/2, Fs/2, N); % Frequency vector centered around zero

for i = 1:3
    F = fft(w(:, i)); % FFT of the signal for each axis
    F = fftshift(F); % Shift zero frequency component to the center
    P(:, i) = abs(F).^2 / N; % Power spectrum
end

% Plot the power spectrum for the Fourier Transform
figure;
semilogy(v, P) % Plots on a log scale on the y-axis
xlim([-1 1]) % Limits the x-axis so as to make the plot more clear
xlabel('Frequency (Hz)', 'Interpreter', 'latex'); % Labels the x-axis
ylabel('Power', 'Interpreter', 'latex'); % Labels the y-axis
title('Power Spectrum', 'Interpreter', 'latex'); % Puts in a title
legend({'\(\omega_x\)', '\(\omega_y\)', '\(\omega_z\)'}, 'Interpreter', 'latex', 'Location', 'best'); % Labels the power spectrum for each of the axis of rotation
grid on;
