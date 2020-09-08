function parabaloid
    clc;
    clear all;
    
    
    max_amplitude = 8.0;
    dx = 100.0;
    dy = 100.0;
    xmax = 8000.0;
    ymax = 8000.0;
    
    scale = 1.0/((xmax/max_amplitude)*sqrt(max_amplitude));
    x = -xmax:dx:xmax;
    y = -ymax:dy:ymax;
    
    [X,Y] = meshgrid(x,y);
    level_floor = zeros(size(X));
    z = -(X.^2/(scale * xmax)^2) - (Y.^2/(scale * ymax)^2) + (max_amplitude);
    
    figure(2)
    hold on
    surf(X,Y,z)
    surf(X,Y,level_floor)
%     zlim([0.0 10.0])
%     view(-37.5,30)
    hold off
    
end