
function [the_x, indvXYRot, Rho]=BestRotation(indvXY, rng) %e.g. rng=[2:15] (Face2~Face15)
    % Return the best angle to rotate the MDS that maximize the cosine
    % similarity between the pairwise angle matrix estimated from the MDS and that estimated from the model  
    % indvXY is Nx2 matrix (X and Y coordinate of each point in MDS)
    % rng is 1D matrix such as [2:15]; To change the model matrix change
    % CosAngRDM.m function
    
            for x0=1:25 %Narrowing down search range in a 15 degree (pi/12) bin
                x=pi/12*(x0-1);
                angM=[cos(x), -sin(x); sin(x), cos(x)];
                newindvXY=indvXY*angM;
                Cs(x0)=CosAngRDM(rng, newindvXY); % return the corr coef between the matrix of pairwise angles estimated from the structural model (rad) and the angles between two points in MDS.
            end
            peak=find(ismember(round(Cs,3), max(round(Cs,3))));
            if numel(peak)>1
                peak=peak(1);
            end
            for a0=1:31 % Narrowing down seach range in a 1 degree bin within the selected 15 degree bin
                x=((peak-2)*pi/12)+((a0-1)*pi/180);
                angM=[cos(x), -sin(x); sin(x), cos(x)];
                newindvXY=indvXY*angM;
                Cs2(a0)=CosAngRDM(rng, newindvXY);
            end
            peak2=find(ismember(round(Cs2,3), max(round(Cs2,3))));
            if numel(peak2)>1
                peak2=peak2(1);
            end
            the_x=((peak-2)*pi/12)+((peak2-1)*pi/180); %rotation angle that maximize Rho
            indvXYRot=indvXY*[cos(the_x), -sin(the_x); sin(the_x), cos(the_x)]; %rotation matrix
            Rho=max(round(Cs2,3)); %maximum correaltion between Euclidean Model and the Test model when rotating the test model in 'the_x' rad.
end