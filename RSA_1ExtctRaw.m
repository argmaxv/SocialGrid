%% RSA ROI
% 1. Extract raw signals from ROIs
% by SPARK 1.Oct.2018

clear;
close all;
clc;

%% Setting

[ProjSet, fs, info, ROI, fname]=Call_default_PS;
raw_path=ProjSet.DATApath;% raw data path
maskpath=ProjSet.ROIpath; %ROI path
[subj, nsubj] = CallSubj_PS(2); %(2) indicates scans acquired from two days are called separately
b0(1)=0; % Stating point of EPI data (cumulative number of EPI scans)
ROIs=ROI.Ana; %ROI list
numROI=numel(ROIs);

%% Main
for roiIdx=1:numROI % number of ROIs
% Select ROI
    ROIselec=roiIdx;
    roitype = ROIs{roiIdx};
    ROIpath=[maskpath, roitype, fs]; % where [roitype]_roi.mat is
    mask=fullfile(ROIpath, [roitype, '_roi.mat']);
    disp(['*** ', roitype, '-', num2str(ROIselec), '/', num2str(numROI), ' ***']);
    
% Extract Raw Timeseries
    for s=1:nsubj

        for r=1:info.Nses
            subpath=[raw_path, subj{s}, fs, info.prefix.Run, num2str(r)];
            disp(['Subject ', num2str(s), ' ', fs,' ', num2str(numel(subj))]);
            clear rawlist Tmean Traw info
            rawlist=dir([subpath, fs, 'wua*.img']);            
            b0(r+1)=size(rawlist,1);

            for b=1:numel(rawlist)
                steps=floor(numel(rawlist)/5);
                data=[subpath, fs, rawlist(b).name];
                [Tmean,Traw, info] = extract_voxel_values(mask,data);
                RAW.Tmean(b+b0(r),1)=Tmean;
                RAW.Traw(b+b0(r),:)=Traw.I.Ya;            
            end % number of Beta.nii, b
            
        end % session

        RAW.id=subj{s};
        RAW.xyz=Traw.I.xyz;
        RAW.mni=Traw.I.mni;
        RAW.roi=info.regions;

        svpath=[raw_path, roitype, fs]; % save path
        if ~exist(svpath, 'dir')
            mkdir(svpath);
        end
        savename=[svpath, subj{s},'_raw.mat'];
        save(savename, 'RAW');
        clear RAW;
        
    end % number of participants, s
end % number of ROIs
