%% RSA ROI
% 2. NoiseNormalization and generate RDM per ROI
% by SPARK 1.Oct.2018

clear;
close all;
clc;

%% Setting

[ProjSet, fs, info, ROI, fname]=Call_default_PS;
raw_path=ProjSet.DATApath;% raw data path
maskpath=ProjSet.ROIpath; %ROI path
Anal1path=ProjSet.ANA1path; %1st level Analysis path to find SPM.mat
Respath=ProjSet.Respath; %save path
[subj, nsubj] = CallSubj_PS;
ROIs=ROI.Ana; %ROI list
nRun=info.Nses; %3
nDay=info.Nday; %2
models={'Mtv_NMatch14', 'Mtv_AllFaces14'}; % Names of 1st level analysis path (See 'design_name' in Multivariate_Level1.m)
%'Mtv_NMatch14' - 14 faces (except for 1 and 16) while matching the number of presentation of each face
%'Mtv_AllFaces14' - All presentations of 14 faces (except for 1 and 16) including the F0,F1, anf F2 events
svoption=1; % 1 if you want to save the results of RDM as mat file. 0 otherwise.

%% Compute RDM
%   With which model (1st level analysis), the patterns of activity would be
%   estimated.
for m=1:numel(models)

    modelname=models{m};
    disp(modelname);

    % Select ROI
    for roiIdx=1:numel(ROIs) 
        ROIselec=roiIdx;
        ROIpath=ROIs{roiIdx};
        
        ROIRDMcollect.crsSesRDM=[];
        ROIRDMcollect.wthSesRDM=[];

        svpath=[Respath, modelname, fs];
        if ~exist(svpath, 'dir')
            mkdir(svpath);
        end    
        disp(['*** the current ROI is ', ROIpath, ' ***']);
        
        for s=1:nsubj 
            clear RDMwth
            
            for dy=1:nDay % 2 days scans
                clear rawdata
                rawdata=load([raw_path, ROIpath, fs, info.prefix.(['day',num2str(dy)]), subj{s}, '_raw.mat']); % where raw signals are extracted from each ROI (RSA_1ExtractRaw)
                Y=rawdata.RAW.Traw;
                cor.xyz=rawdata.RAW.xyz;
                cor.mni=rawdata.RAW.mni;
                clear SPM spmpath design
                spmpath=[Anal1path, modelname, fs, info.prefix.(['day',num2str(dy)]), subj{s}, fs, 'SPM.mat'];
                design=load(spmpath); %using the SPM.mat generated by 1st level analysis
                SPM=design.SPM;
                [u_hat,resMS,Sw_hat,beta_hat,shrinkage,trRR]=rsa.spm.noiseNormalizeBeta(Y,SPM); %noise normalization
                nBeta=numel(SPM.Sess(1).Fc); % number of Beta in each block
                nEtc=6; %numMotionregressor;
                
    % 1. Compute the within block RDM (for MDS)
                clear Uhat0 Uhat rdm_eachses rdm_eachses4Mean
                for nR=1:nRun
                    stp = 1 + (nR-1)*(nBeta+nEtc);
                    enp = stp + nBeta -1;
                    Uhat0{nR}=u_hat(stp:enp,:);
                    for zi=1:size(Uhat0{nR},1)
                        Uhat{nR}(zi,:)=zscore(Uhat0{nR}(zi,:));
                    end
                    Uhats(nR,:,:)=Uhat{nR};
                    if sum(sum(isnan(Uhat{nR})))~=0
                        error(['There is NaN voxels in run ', num2str(nR)]);
                    end
                    rdm_eachses{nR}=squareform(pdist(Uhat{nR},'Euclidean')); % Within blcok RDM 
                    rdm_eachses4Mean(nR,:,:)= rdm_eachses{nR};
                end
                RDMwth.(ROIpath).RDM(dy,:,:)=squeeze(mean(rdm_eachses4Mean));
                
    % 2. Compute the RDM by comparing different blocks (cross runs; for RSA)
                clear Othhat crsSesRDM crsSesRDM4mean
                for nR=1:nRun
                    clear otherUhat
                    runs=[1:nRun];
                    otherrun=runs(~ismember(runs,nR));
                    for ori=1:numel(otherrun)
                        otherUhat(ori,:,:)=Uhat{ori};
                    end
                        Othhat{nR}=squeeze(mean(otherUhat)); % Mean Beta of others blocks
                        crsSesRDM{nR}=pdist2(Uhat{nR},Othhat{nR},'Euclidean');
                        crsSesRDM4mean(nR,:,:)=crsSesRDM{nR};
                end                
                RDM.(ROIpath).(['RDMx', num2str(dy)]){s}=squeeze(mean(crsSesRDM4mean));
                RDM.(ROIpath).(['RDMx', num2str(dy),'_col'])(:,:,s)=squeeze(mean(crsSesRDM4mean));
            end %for days
            % Average within block RDMs across days per subject (symmetrical; diag=0)
            MDS.(ROIpath).RDM(s,:,:)=squeeze(mean(RDMwth.(ROIpath).RDM));
        end %for Subject
                
        if svoption
            save(fullfile(svpath, fname.RSAroi),'RDM');
            save(fullfile(svpath, fname.MDSroi),'MDS');
        end
    end %for ROIs
end %for Models