function [on, cat] = isRGridN(dgr, n)
% R=radian
% it returns on/off grid as the first value
% it returns the category it belongs among 12 categories with 30 degree bins (if n=6) 
    angle=pi/n;
    totalcat=2*n-1;
    for t=1:length(dgr)
        
        if dgr(t)<0
            dgr(t)=2*pi+dgr(t);
        end
        
        ans=round(dgr(t)/angle);
        if ans==totalcat+1
            ans=0;
        end
        cat(t)=ans;
        on(t)=mod(cat(t),2)==0; % Category 0 2 4 6 ... is on-grid 
    end
end