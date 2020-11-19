
function [ProjSet, fs, info, ROI, fname]=Call_default_PS
    
    % Defult setup for the 'Partner selection' task analysis
    % by APark

%% Basic infomation
    fs=filesep;
    info.rngset='default'; % Random seed
    info.Nday=2; % Number of days
    info.Nses=3; % Number of blocks in a day
    info.Nface=16; % number of entities (faces in the social hierarchy) 
    info.Nperm=1000; % Number of permutation
    info.tunit=1000; % Unit of reaction time 1s=1000ms
    info.periodicity=6; %6 fold periodicity
    info.prefix.day1='P';
    info.prefix.day2='Q';
    info.prefix.Run='R'; %img/hdr files are in the folder, R1, R2, and R3 according to in which block there were acquired.
    info.PhiModel='Grid6_F01_F02_F12_FDec_5s';
%% ROIs
% ROIs should have the same name (in ROIpath)
    ROI.Ana={'HCr', 'HCl', 'ECr', 'ECl', 'M1r', 'M1l'}; % Anatomically defined ROIs
    ROI.GP={'mPFC_GP', 'mPFC_GP', 'mPFC_GP'}; % Based on the univariate analysis decision value (GP)
    ROI.Grid={'EC_Grid', 'mPFC_Grid'}; % Based on the F-test for hexadirectional grid-like codes
    
%%  Base path - update them according to your system
    projpath='/home/ravic/Pr4.PF/';
    ProjSet.basepath=[projpath, 'Imaging'];
    
    %The following programs/toolboxes are needed
    ProjSet.spmdir='/usr/local/MATLAB/spm12';
    ProjSet.rsatoolbox=[projpath, 'Programs/rsatoolbox-develop'];
    ProjSet.circstat=[projpath, 'Toolbox/circstat-matlab-master'];
    % diedrichsenlab_util
    % marsbar-0.44
    % surfing-master
    
    addpath(ProjSet.basepath);
    addpath(ProjSet.spmdir);
    addpath(genpath(ProjSet.rsatoolbox));
    addpath([ProjSet.basepath, fs, 'Batch']);
    
%%  Analyses path
    ProjSet.func=[ProjSet.basepath, fs, 'Batch', fs, 'func']; % Where custom functions are
	if ~exist('func','dir')
        addpath(ProjSet.func);
    end
    %behavior
    ProjSet.ONSETpath    = [ProjSet.basepath, fs, 'Onset', fs]; % path to save behavioral data (mat)
    ProjSet.Logpath         =[ProjSet.ONSETpath, 'fMRI_Bhv_logtxt', fs]; % path where behavioral data are (txt log files)
    ProjSet.Metapath       =[ProjSet.ONSETpath, 'fMRI_Bhv_mat', fs]; % path to save behavioral matadata
    %imaging
    ProjSet.DATApath      =[ProjSet.basepath,  fs, 'Data', fs];    % path where raw neuroimaging data are
    ProjSet.MODELpath   =[ProjSet.basepath,  fs, 'Batch', fs, 'Model', fs]; % where model and contrast are defined for GLMs
    ProjSet.ANApath        =[ProjSet.basepath,  fs, 'Analysis', fs]; % base path for analysis
    ProjSet.ANA1path      =[ProjSet.ANApath, 'Analysis Lv1', fs]; %lv1 analysis
    ProjSet.ANA2path      =[ProjSet.ANApath,  'Analysis Lv2', fs]; %lv2 analysis (univariate analyses)
    ProjSet.Phipath          =[ProjSet.ANApath, 'Analysis Phi', fs]; %path to save grid analysis
    %ProjSet.PhiInfopath    =[ProjSet.Phipath, 'CurrentModel', fs]; %path where phi inomation is
    ProjSet.PhiInfopath     =[ProjSet.Phipath, 'GeneratedModel', fs]; %path where phi inomation is
    ProjSet.ROIpath         =[ProjSet.ANApath, 'Analysis ROI', fs];  %ROI path
    ProjSet.Maskpath      =[ProjSet.ROIpath, 'Mask100', fs]; 
    ProjSet.Respath        =[ProjSet.ANApath, 'Analysis RSA', fs]; % Save folder for multivariate analyses
    
    addpath('/home/ravic/Pr4.PF/Toolbox/circstat-matlab-master');
    
%%  Save file names
    fname.savename='sl_RSA_tauA.mat';
    fname.savefolder='TauMap_Euc';
    fname.maskname='Mask.img';
    fname.searchlightname='searchlight_100.mat';
    
end
