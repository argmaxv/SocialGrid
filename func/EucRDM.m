function [Model_rdm, Model_rdm_perm]=EucRDM(rng, Nperm)% [2:15]';
    
    rand('twister',0);
    q=[1:4];
    rng=rng';
    [X,Y]=meshgrid(q,q);
    X=X(rng);
    Y=Y(rng);
    I=[1:length(X)];
    [i,j]=meshgrid(I,I);
    Model_rdm=sqrt((X(i)-X(j)).^2+(Y(i)-Y(j)).^2);
   
    if nargin>1
        Xperm=randi([1,4],length(X),Nperm);
        Yperm=randi([1,4],length(Y),Nperm);

        for n=1:Nperm
            clear X0 Y0
            X0=Xperm(:,n);
            Y0=Yperm(:,n);
            I=[1:length(X)];
            [i,j]=meshgrid(I,I);
            Model_rdm_perm{n,1}=sqrt((X0(i)-X0(j)).^2+(Y0(i)-Y0(j)).^2);
        end
    end
end
    