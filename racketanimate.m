close all;
% Define the semi-major and semi-minor axes lengths
a = 3; % Semi-major axis length
b = 2; % Semi-minor axis length
c = 1; % Third semi-axis length (for 3D ellipse)
% Create a 3D ellipse
theta = linspace(0, 2*pi, 100);
z = c * ones(size(theta)); % Constant z-coordinate for the ellipse
x = a * cos(theta);
y = b * sin(theta);

s = plot3(x,y,z);

for i = 0:360
    rotate(s,[1 1 1],5);
    drawnow;
    pause(0.1);
    axis([-3 3 -3 3 -3 3]);
end
% this spins. next step, add a marker for handle (line?)
% incorporate into main.m using updated rotation matrix
%%
close all;
% Define the semi-major and semi-minor axes lengths
a = 3; % Semi-major axis length
b = 2; % Semi-minor axis length
c = 1; % Third semi-axis length (for 3D ellipse)
% Create a 3D ellipse
theta = linspace(0, 2*pi, 100);
z = c * ones(size(theta)); % Constant z-coordinate for the ellipse
x = a * cos(theta);
y = b * sin(theta);

% Plot the ellipse
s = plot3(x,y,z, 'b');
hold on;
% Plot a marker at the start of the ellipse
plot3(x(1),y(1),z(1),'r*');
hold off;

% Group the ellipse and the marker for rotation
h = hgtransform;
set(s, 'Parent', h);
% Now rotate the group
for i = 0:360
    rotate(h,[1 1 1],5);
    drawnow;
    pause(0.1);
    axis([-3 3 -3 3 -3 3]);
end
% this spins. next step, add a marker for handle (line?)
% incorporate into main.m using updated rotation matrix