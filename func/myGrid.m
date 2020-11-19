function [rad, euc] = myGrid(x, y)
% [x,y]=[X1,Y1]-[X0,Y0]
% This function returns radian degree and euclidian distance
    for t=1:length(x)
        rad(t)=(atan2(y(t),x(t)));
        euc(t)=sqrt(x(t)^2+y(t)^2);
    end
end