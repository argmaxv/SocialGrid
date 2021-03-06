%% RSA searchlight 
% 1. Extracting Vols from spmT and generate RDM per searchlight
% by Alex Park 10.Oct.2018

close all
clear
clc

%% turn on the parallel computing
pflag = gcp('nocreate');
if isempty(pflag)
    poolsize=0;
else
    poolsize=pflag.NumWorkers;
end

%% Path setup
[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[Subj, nSubj]=CallSubj_PS(2); %Define subject names sepaarately if it was acquired from differnt days
Respath=ProjSet.Respath;
BrRDM={'Mtv_FaceNMatch14', 'Mtv_AllFaces14'}; % Select the activity patterns generated by which model to generate RDM
setups.Lv1path=ProjSet.ANA1path;
setups.Subj=Subj;
setups.NSub=nSubj;
setups.NSes=info.Nses; % Number of sessions
Vmask=spm_vol(fullfile(ProjSet.Maskpath, fname.maskname)); %Mask [1, inside of the brain; 0, outside of the brain]
Vmask.data=spm_read_vols(Vmask);
mask=nan(size(Vmask.data));
mask(Vmask.data==1)=1;
setups.mask=mask;
setups.searchlight=fullfile(ProjSet.Maskpath, fname.searchlightname); % Searchlight defined from RSA_sL0_GenSearchlight.m
setups.RDMsearchligh=fname.RDMsearchligh;

%% Main
for Bri=1:numel(BrRDM)
    disp(BrRDM{Bri});
    savepath=[Respath, BrRDM{Bri}, fs, 'Vols']; %save path
    if ~exist(savepath,'dir')
        mkdir(savepath);
    end
    setups.savepath=savepath;
    setups.model=BrRDM{Bri};
    tic;
    subjectPar(setups);
    toc;
end 
delete(gcp('nocreate'));
    
function subjectPar(setups)
    for Sbi=1:setups.NSub
        disp(setups.Subj{Sbi});
        RDMpath=[setups.Lv1path, setups.model, filesep, setups.Subj{Sbi}];
        spmTlist=dir([RDMpath, filesep, 'spmT*.nii']);
        spmTname={spmTlist.name};
        Tmat = extractT(RDMpath, spmTname, setups); % Extract beta activity from t map
        RDMs = makeRDM(Tmat, setups);   % Generate RDM. the dissimilarity is estimated by Euclidean distance between activity patterns (cross-session) 
        save([setups.savepath, filesep, setups.Subj{Sbi}, setups.RDMsearchligh],'RDMs');
        clear RDMs    
    end %Sbj number
end

function Tmat = extractT(RDMpath, spmTname, setups)
% Extracting beta signals from a T map and whole brain zscore transform
    parfor spmTi=1:numel(spmTname)
        spmTpath=[RDMpath, filesep, spmTname{spmTi}];
        data = spm_read_vols(spm_vol(deblank(spmTpath)));
        maskeddata=data.*mask;
	    zdata=(maskeddata-nanmean(maskeddata(:))) / nanstd(maskeddata(:));
	    %zdata=nanzscore(data.*mask);
        Tmat(spmTi,:) = zdata(:);
    end
end

function [RDMs] = makeRDM(Tmat, setups)
% Segmenting beta series per serchlight according to their coordinates
    load(setups.searchlight); %Ll, which include the coordinate of 100 voxels of each searchlight
    nROIs=size(LI,1);    
    parfor LIi=1:nROIs
        SingleSubVol=Tmat(:,LI{LIi});
        RDMs(LIi).RDM =genRDM(SingleSubVol,setups.NSes);
    end
end

function [RDMs]=genRDM(SingleSubVol,nses)
        nCondition=size(SingleSubVol,1)/nses; %RDM size
        % activity patterns of the target session
        for s=1:nses
            SingleSubVol0=SingleSubVol(nCondition*(s-1)+1:nCondition*s,:);
            SingleSubVol_ses{s}=ZinROI(nCondition, SingleSubVol0);
        end
        % mean activity patterns of the other sessions
        sesmat=[1:nses];
        for s=1:nses
            altses=sesmat(~ismember(sesmat,s));
            for alts=1:numel(altses)
                SingleSubVol_alt_temp(alts,:,:)=SingleSubVol_ses{altses(alts)};
            end
            SingleSubVol_alt{s}=squeeze(nanmean(SingleSubVol_alt_temp));
        % Compute RDM from the Euclideandistance between activity patterns
        % (from the target to other sessions)
            RDM(s,:,:)=pdist2(SingleSubVol_ses{s}, SingleSubVol_alt{s} , 'euclidean');
        end
        RDMs =squeeze(nanmean(RDM));
end

function [SingleSubVol]=ZinROI(nCondition, SingleSubVol0) % zscore
    for i=1:nCondition
        SingleSubVol(i,:)=nanzscore(SingleSubVol0(i,:));
        %SingleSubVol(i,:)=(SingleSubVol0(i,:) - nanmean(SingleSubVol0(i,:)) ) / nanstd(SingleSubVol0(i,:));
    end
end
