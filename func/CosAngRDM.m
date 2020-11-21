function [Cs, Cs_perm, modelRDM, testRDM]=CosAngRDM(rng, XY, Nperm)
% Retun the correlation of pairwise angles between Euclidean_model
% (modelRDM) and test_model (testRDM, e.g. MDS)
% Inputs: 
%e.g. rng=[2:15] Face2~Face15
%XY = x and y coordinate (e.g. MDS) 2xnumel(rng) matrix
%NPerm, num permutation

    rand('twister',0);
    q=[1:4];
    rng=rng';
    [X,Y]=meshgrid(q,q);
    X=X(rng);
    Y=Y(rng);
    I=[1:length(X)];
    [i,j]=meshgrid(I,I);
    %Model_rdm=sqrt((X(i)-X(j)).^2+(Y(i)-Y(j)).^2);
    modelRDM=cos(atan2( (Y(j)-Y(i)), (X(j)-X(i)) ));
    
     Y2=XY(:,2);
     X2=XY(:,1);
     testRDM=cos(atan2( (Y2(j)-Y2(i)), (X2(j)-X2(i)) ));    
    Cs = corr(modelRDM(:),testRDM(:));

    if nargin>2
        XYrng=(normalize(XY,'range')* (max(rng)-min(rng)) )+min(rng);
        Y0=XYrng(:,2);
        X0=XYrng(:,1);
        Xperm=(max(rng)-min(rng)).*rand(length(X), Nperm)+min(rng);
        Yperm=(max(rng)-min(rng)).*rand(length(Y), Nperm)+min(rng);

        for n=1:Nperm
            clear X0 Y0 testRDM_perm
            X0=Xperm(:,n);
            Y0=Yperm(:,n);
            I=[1:length(X)];
            [i,j]=meshgrid(I,I);
            testRDM_perm=cos(atan2( (Y0(j)-Y0(i)), (X0(j)-X0(i)) ));
            Cs_perm(n) = corr(testRDM_perm(:),testRDM(:)); 
        end        
    end    
end
    