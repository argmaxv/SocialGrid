function [a] =faceplace(face,type,sub)
% Get the position [x,y] of face in differnt maps.
% [a] =faceplace(face,type,sub)
% i.e. faceplace(3,'p',2)
    for lnth=1:numel(face)
        if nargin==1 || (nargin==3 && strcmp(type,'e')) %euc
            pos=[1,1;2,1;3,1;4,1;1,2;2,2;3,2;4,2;1,3;2,3;3,3;4,3;1,4;2,4;3,4;4,4];
            a(lnth,1)=pos(face(lnth),1);
            a(lnth,2)=pos(face(lnth),2);
        elseif nargin>1
            
            if strcmp(type,'p') %placement
                clear pos
                pos=load('/home/ravic/Pr4.PF/Placement/Placement_new/place_coor25.mat');
                a(lnth,1)=pos.place.x(sub,face(lnth));
                a(lnth,2)=pos.place.y(sub,face(lnth));
        
            elseif strcmp(type,'s') %sigmoid %x=[-3:2:3]; y=1+(3./(1+exp(-x)))
                clear pos
                pos=[1,1;2,1;3,1;4,1;1,2;2,2;3,2;4,2;1,3;2,3;3,3;4,3;1,4;2,4;3,4;4,4];
                pos=pos*2-5;
                possigmoid=1+(3./(1+exp(-pos)));
                a(lnth,1)=possigmoid(face(lnth),1);
                a(lnth,2)=possigmoid(face(lnth),2);

            elseif strcmp(type,'t') %Tangent 
                clear pos
                pos=[1,1;2,1;3,1;4,1;1,2;2,2;3,2;4,2;1,3;2,3;3,3;4,3;1,4;2,4;3,4;4,4];
                pos=pos*2-5;
                posTan=tan(pos);
                a(lnth,1)=posTan(face(lnth),1);
                a(lnth,2)=posTan(face(lnth),2);     
            
            elseif strcmp(type,'m') %MDS position
                clear pos pos0
                %pos=load('/home/ravic/Pr4.PF/Placement/Placement_new/mds_coor25.mat');
                %pos0=load('/home/ravic/Pr4.PF/Placement/Placement_new/mds_rot_pos25.mat'); %rotation and x<->y
                %pos0=load('/home/ravic/Pr4.PF/Placement/Placement_new/mds_Rotation16_coordinate.mat'); %16 all rotation - /home/ravic/Pr4.PF/Placement/MDS/Mtv_DMa16_test/rng_1_16
                pos0=callMDS16;
                pos=pos0.mdspos.HC; %=pos0.mdspos.EC
                a(lnth,1)=pos.place.x(sub,face(lnth));
                a(lnth,2)=pos.place.y(sub,face(lnth));
                
            else
                error('error in Face->Value in faceplace.m. Check your code');
            end
            
        end%if1
    end%for
end%func

