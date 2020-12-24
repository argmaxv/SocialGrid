%% RSA searchlight 
% 3. Map the ranks correation to the center of each searchlight and
% smoothing
% by Alex Park 10.Oct.2018

%	read KendallT (T.(RDM) has sub x voxels matrix)
%	smooth individual tau map
%	save the smoothed tau map to 1d mat (sKendallT.mat) and nii file ('s')
%   write to nii images (MaskMap size)

close all
clear
clc
%% Turn on the parallel computing
    pflag = gcp('nocreate');
    if isempty(pflag)
        poolsize=0;
    else
        poolsize=pflag.NumWorkers;
    end
%% Path setup
[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[Subject, NSubject]=CallSubj_PS; 

% Define the Brain RDM path
ModelNames={'Mtv_FaceNMatch14', 'Mtv_AllFaces14'};
BrDataPath=ProjSet.Respath;

% Mask info to make 1D array to 3D map
Mask = fullfile(ProjSet.Maskpath, fname.maskname);
struct=spm_vol(deblank(Mask));
MaskMap =spm_read_vols(struct);
vsize=struct.dim; 
load(fullfile(ProjSet.Maskpath, fname.searchlightname)); %Ll, the coordinates of each searchlight
prefix=info.prefix.smooth;
FWHMmm=info.FWHMmm;

%% Main
for mm=1:numel(ModelNames)
    ModelName=ModelNames{mm};
    datapath=[BrDataPath, ModelName, fs];
    savepath=[BrDataPath, ModelName, fs, 'sL', fs];
    mtype = fname.RSAsl;
    load(fullfile(savepath, mtype)); %T
    RDMnames=fieldnames(T); % Subjects
    svfolder='TauMap_Euc_T';
    nRDM=size(RDMnames,1); % Num Subjects
    for iRDM=1:nRDM
        sSimilarity=[];
        disp(RDMnames{iRDM});
        mappath=fullfile([savepath, svfolder], RDMnames{iRDM});
        if ~exist(mappath)
            mkdir(mappath);
        end
        parfor sub=1:size(T.(RDMnames{iRDM}),1)
            thearray = T.(RDMnames{iRDM})(sub,:);
            filename = [ModelName, '_', num2str(sub),'_', RDMnames{iRDM}];
            [~, fname] = array2nii(thearray, voxel, MaskMap, struct, mappath, filename); %save unsmoothed nii
            smoothingnii(FWHMmm, prefix, fname, 1); %save smoothed nii withprefix 's'
        end
    end
    disp('Done!');
end
