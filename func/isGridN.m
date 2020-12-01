function [on cat] = isGridN(x, y, n)
% [x,y]=[X1,Y1]-[X0,Y0]
% compute the cos(theta) and radien degree
% it returns on/off grid as the first value
% it returns the category it belongs among 12 30 degree bins  
    angle=180/n;
    totalcat=2*n-1;
    
    for t=1:length(x)
        %dgc(t)=cos(atan2(y(t),x(t)));
        dgr(t)=atan2(y(t),x(t))*180/pi;
        
        if dgr(t)<0
            dgr(t)=360+dgr(t);
        end
        ans=round(dgr(t)/angle);
        if ans==totalcat+1
            ans=0;
        end
        cat(t)=ans;
        on(t)=mod(cat(t),2)==0; 
    end
end