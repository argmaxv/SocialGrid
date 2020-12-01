%% EachDec_Beta_collect.m
    % Get the mean beta signals from the ROI during presentation of each pair (F0-F1 and F0-F2)
    % Save the mean beta (num. pairs x num. subjects), the beta values of each participants are
    % z-scored within each block.

clear; close all; clc;
[ProjSet, fs, info, ROI, fname]=Call_default_PS;
%%%%%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%
% option1 raw activity
model_name='EachDec';                       % 1st lev GLM model in which the Beta files are generated
ROIs=[ROI.Grid, ROI.GridAna];                % 'EC_Grid', 'mPFC_Grid', 'ECr', 'mPFC' (In which the mean betas are extracted)
PM_interest='';                                       % raw activity, no parametric modulations
% % option2, alternatively (testing GP effects)
% model_name='EachDec_GP_6folds'; %model_name='EachDec_GP';             % 1st lev GLM model in which the Beta files are generated
% ROIs=ROI.GP;                                      % 'mPFC_GP', 'mPFC_GP', 'mPFC_GP'
% PM_interest='GP';

% Set which Phi to apply
phiset.ROI = ROI.Grid{1};                       % which grid angle will extract (e.g. EC_Grid), note that this is not the ROI where extract the beta signals.
phiset.type = 'PhixB';                              % cross validation of Phi which is estimated differnt blocks (across blocks); PhixD - across days 
svoption=1;                                             % 1 to save the results 0 otherwise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sub0, subn] = CallSubj_PS;
BrainDataPath=[ProjSet.ANA1path, model_name, fs];
ROIpath=ProjSet.ROIpath;
%% Get the labels of each Beta from SPM.m (1st lv GLM specification)
[BetaInfoPM, BetaPMNum, ~]=Get_Beta_Labels(model_name, PM_interest);
%Returns 'BetaInfoPM', 'BetaPMNum'

%% Categorizing Each subject's Theta-Phi
    % See GridCategory.m
%[AltPeriodicity, ~]=Altanative_Periodicity(phiset); % See GridCategory.m 

%% Turn on the parallel computing
    pflag = gcp('nocreate');
    if isempty(pflag)
        poolsize=0;
    else
        poolsize=pflag.NumWorkers;
    end
    
%% Extract Betas

for ridx=1:numel(ROIs)
    clear curroi nonpmB zB
    curroi=fullfile([ROIpath, fs, ROIs{ridx}], [ROIs{ridx}, '_roi.mat']);
    fprintf([ROIs{ridx}, '\n']);
    
    parfor sub=1:subn
        fprintf([sub0{sub}, '\n']);
        cursubdir=[BrainDataPath, fs, sub0{sub}]; % Each participants folder in 1st lv analysis where beta files are
        nonpmB(:,sub)=getBetas(cursubdir, curroi, BetaPMNum); %BetaPMNum is the index of Beta file to extract which is generated from Get_BetaInfo.m
    end
    beta_block=[BetaInfoPM(:).Block]';
    beta_names={BetaInfoPM(:).EventType}';
    parfor sub=1:numel(sub0)
        zB(:,sub)=ZBetaperBlock(nonpmB(:,sub), beta_block); %Beta files are z scored within each block
    end
    if svoption
        ExtractedBetaPath=[ProjSet.PhiInfopath, phiset.ROI, fs, 'zB', fs]; % where to save the results
        if ~exist(ExtractedBetaPath, 'dir')
            mkdir(ExtractedBetaPath);
        end
        save(fullfile(ExtractedBetaPath, [ROIs{ridx}, '.mat']), 'beta_block', 'beta_names', 'zB');
        %save(fullfile(ExtractedBetaPath, ['GridCat_', ROIs{ridx}, '_',
        %phiset.type, '.mat']), 'AltPeriodicity'); %GridCat_EC_Grid_PhixB
        %(This part is done in GridCategory.m)
    end
end
