function [E_rdm, D1_rdm, D2_rdm, E_rdm_perm, D1_rdm_perm, D2_rdm_perm]=EucRDM_factorize(range, Nperm)% [2:15]';E_rdm, 
    
    %rand('twister',0);
    rng('default');
    q=[1:4];
    range=range';
    [X,Y]=meshgrid(q,q);
    X=X(range);
    Y=Y(range);
    I=[1:length(X)];
    [i,j]=meshgrid(I,I);
    E_rdm=sqrt((X(i)-X(j)).^2+(Y(i)-Y(j)).^2);
    %GP=max(X(i), X(j)).*max(Y(i),Y(j));
    D1_rdm=abs(X(i)-X(j));
    D2_rdm=abs(Y(i)-Y(j));
   
    if nargin>1
        Xperm=randi([1,4],length(X),Nperm);
        Yperm=randi([1,4],length(Y),Nperm);

        for n=1:Nperm
            clear X0 Y0
            X0=Xperm(:,n);
            Y0=Yperm(:,n);
            I=[1:length(X)];
            [i,j]=meshgrid(I,I);
            E_rdm_perm{n,1}=sqrt((X0(i)-X0(j)).^2+(Y0(i)-Y0(j)).^2);
            D1_rdm_perm{n,1}=abs(X0(i)-X0(j));
            D2_rdm_perm{n,1}=abs(Y0(i)-Y0(j));
        end
    end
end
    