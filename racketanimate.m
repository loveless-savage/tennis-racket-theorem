%% Animation Station
dt = 0.1;
te = 0:dt:tmax; % evenly spaced time array
we = interp1(t,w(:,1),te,'spline'); % fix when we know which column of w is theta

racket = images.roi.Ellipse('SemiAxes',[10 20],'RotationAngle',0);
for istep=1:length(te)
    
    axis square
    pause(dt-0.01)
end

plot(racket)