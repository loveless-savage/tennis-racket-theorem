% demonstration that the rotation matrices work correctly:
% multiply each unit vector by the rotation matrix of each frame, then
% watch the endpoints trace out curves on the unit sphere
figure;
for n=1:N
    clf;
    plot3(shiftdim(rot(1,1,1:n)),shiftdim(rot(2,1,1:n)),shiftdim(rot(3,1,1:n)));
    hold on; grid on;
    plot3(shiftdim(rot(1,2,1:n)),shiftdim(rot(2,2,1:n)),shiftdim(rot(3,2,1:n)));
    plot3(shiftdim(rot(1,3,1:n)),shiftdim(rot(2,3,1:n)),shiftdim(rot(3,3,1:n)));
    axis([-1 1 -1 1 -1 1]);
    title('Rotation Integration Test');
    subtitle('Tracking each unit vector over time');
    xlabel('x');ylabel('y');zlabel('z');
    legend('x','y','z',Location='northeast');
    pause(0.01);
end