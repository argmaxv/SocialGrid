% 2D MDS for visualization

clear 
close all;
clc

%% Input model
[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[subj, subn] = CallSubj_PS;
model={'Mtv_NMatch14', 'Mtv_AllFaces14'}; % See model names in RSA_2NoiseCov.m
%, 'Mtv_AllFace16'}; 
ROIs=ROI.Ana;
svoption=0; %1 to save the results 0 otherwise
NPerm=info.Nperm; % number of permutation - test how mean group MDS capture structure better compared to the random 14 points in the same areas.

%% main
for m=1:numel(model)
    clear ROIRDM.ComSesBi;
    fprintf('\n %d. %s ', m, model{m});
    sz=model{m}(end-1:end);
    Face_rng=1:str2num(sz);
    if str2num(sz)==14 %model={'Mtv_NMatch14', 'Mtv_AllFaces14'}
        RDMModel=EucRDM(Face_rng+1); %Face 2:15 were included in MDS
    elseif str2num(sz)==16 %model={'Mtv_AllFace16'};
        RDMModel=EucRDM(Face_rng); %Face 1:16 were included in MDS
    end    
    
    clear MDS
    MDSpath=[ProjSet.Respath, fs, model{m}, fs];
    if svoption
        svpath=[MDSpath, 'MDS'];
        savefilename=['MDS_Results_', model{m}];
        if ~exist(svpath, 'dir')
            mkdir(svpath)
        end
    end
    RDM0path=fullfile(MDSpath, fname.MDSroi);            
    load(RDM0path); %MDS generated from RSA_2NoiseCov.m
    
    for i=1:2:4 %2: combine right and left MDS; 4: numROIs - HCr,HCl,ECr, and ECl;
        ROIName=ROIs{i}(1:end-1);
        fprintf('\n ROI=%s \n', ROIName);
	%% individual MDS
        for s=1:subn
                fprintf('%s ', '.');
                clear pairwsdist 
                pairwsdist_r=squeeze(MDS.([ROIName,'r']).RDM(s,:,:));
                pairwsdist_l=squeeze(MDS.([ROIName,'l']).RDM(s,:,:));
                pairwsdist=(pairwsdist_r+pairwsdist_l)/2;
                ROIRDM.(ROIName).RDMBi(s,:,:)=pairwsdist;   
                [mds_xy0, ~]=mdscale(pairwsdist, 2); %2-D
                dim1(:,s)=mds_xy0(:,1);
                dim2(:,s)=mds_xy0(:,2);

    %% Find rotate matrix for each individual MDS that maximize corr(Rho) with that of the 4x4 structural model (Euclidean Model)
     %  the rotated MDS is only used to find the brain activity encoding
     %  hex. grid-like codes modulated by the inferred trajectories defined
     %  on MDS coordinate.
                clear indvXY;
                indvXY=[dim1(:,s), dim2(:,s)];
                [the_x, indvXYRot, RhoMax]=BestRotation(indvXY, Face_rng);            
                indv.(ROIName).dim(:,:,s)=indvXYRot; % X,Y coordinate in a rotated matrix (XY x RotationMatrix)
                indv.(ROIName).rad(s)=the_x; % Roation angle
                indv.(ROIName).Rho(s)=RhoMax; %Max Pearson correation between the pairwise cos angles of Euclidean Model and Test(MDS) model  
        end %for subj s
        
	%% MDS figure based on a group mean
        fprintf('\n ');
        mdsfig=figure;
        ROIRDM.(ROIName).ComBi= squeeze(mean(ROIRDM.(ROIName).RDMBi));
        %RDM0mean=squeeze(mean(ROIRDM.(ROIName).ComBi,1)); % generating a group level dissimilarity matrix based on mean of subjects computed in each cell
        [mds_XY0, ~]=mdscale(ROIRDM.(ROIName).ComBi,2);
        mds_XY=mds_XY0; %*RotM; %no rotation
        c = linspace(1,10,length(mds_XY)); %point colors
        mdsfig=scatter(mds_XY(:,1), mds_XY(:,2), length(Face_rng), c, 'filled');
        hold on;
        for idx=1:length(mds_XY) % Putting the text labels
            text(mds_XY(idx,1),mds_XY(idx,2),num2str(Face_rng(idx)));
        end
        title([ROIName ' MDS ' model{m}]);
        mdsXY.(ROIName).XY=mds_XY0;

        if svoption
            saveas(mdsfig, fullfile(svpath, ['GroupMDS_', ROIName]), 'epsc');
            saveas(mdsfig, fullfile(svpath, ['GroupMDS_', ROIName,'.fig']) );
            close(gcf);
        end

    %% Permutation test
        clear mdsRDM mdsRDMPerm cs_mds cs_perm
        [mdsRDM, mdsRDMPerm]=mdsRDM(mds_XY0, Face_rng, NPerm); % Make a pairwise distance matrix based on MDS points.
        for p=1:NPerm
            if mod(p,NPerm/20)==0 fprintf('%s ', '.'); end
            Rho.(ROIName).Perm(p)=corr(RDMModel(:), mdsRDMPerm{p}(:), 'type', 'Spearman'); %Baseline distribution of correation coeffieicnt. 
        end
        h=0; rdbu=RdbuMap;
        h=1+h; rdmfig(h)=figure; imagesc(RDMModel); colormap(rdbu); title([model{m}, ' EUC dist ', ROIs{i}]); colorbar; 
        h=1+h; rdmfig(h)=figure; imagesc(mdsRDM); colormap(rdbu); title([model{m}, ' MDS dist ', ROIs{i}]); colorbar; 

        Rho.(ROIName).mds=corr(RDMModel(:), mdsRDM(:), 'type', 'Spearman');
        Rho.(ROIName).Perm(NPerm+1)=Rho.(ROIName).mds;
        Rho.(ROIName).pVal=sum(Rho.(ROIName).Perm>Rho.(ROIName).mds)/(NPerm+1);

        [rotang, indvXYRot0]=BestRotation(mds_XY0, Face_rng); 
        [cs_mds, cs_perm, AgRDM, AgMDS]=CosAngRDM(Face_rng, indvXYRot0, NPerm);
        h=1+h; rdmfig(h)=figure; imagesc(AgRDM); colormap(rdbu); title([model{m}, ' Euc Cos ', ROIs{i}]); colorbar; 
        h=1+h; rdmfig(h)=figure; imagesc(AgMDS); colormap(rdbu); title([model{m}, ' MDS Cos ', ROIs{i}]); colorbar; 

        CorrAng.(ROIName).mds=cs_mds;
        CorrAng.(ROIName).Perm=[cs_perm, cs_mds];
        CorrAng.(ROIName).pVal=sum(CorrAng.(ROIName).Perm>CorrAng.(ROIName).mds)/(NPerm+1);
        if svoption
            savefig(rdmfig, fullfile(svpath, [savefilename,'_', (ROIName), '.fig']) );
            close(rdmfig);
        end
    end % for ROI r 
    if svoption
        if ~exist(svpath, 'dir')
            mkdir(svpath)
        end
        save(fullfile(svpath, [savefilename, '.mat']), 'MDS', 'Rho', 'CorrAng', 'mdsXY', 'indv');
    end
    fprintf('\n ');
end %for model
