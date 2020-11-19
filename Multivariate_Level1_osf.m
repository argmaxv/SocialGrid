% 1st level analysis for multivariate analysis.
%
% by SPARK 12/28/2017

clear
close all;
clc;

%% Path setup

[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[subj, nsubj] = CallSubj_PS(2); %(2) indicates that the scan aqcuired in differnt days will be modeled separately in this GLM (see CallSubj_PS.m) 
nblocks      = info.Nses;
base_path =ProjSet.basepath; 
data_path =ProjSet.DATApath;
an_path =[ProjSet.ANA1path];            % Where results of lv1 analysis are
ons_dir=ProjSet.ONSETpath;               % Behavioral data
addpath (ProjSet.MODELpath);
motion_reg=6;                                    % number of motion regressoes (default is 6);
flag=0;                                                 % for sanity check (if flag=1 then error message)

%%%%%%%%%%%%%%%%%%% Input the following information %%%%%%%%%%%%%%%%%%%%%%%%%

% 1. model name - choose the model according to which data should be included into the model 
design_name = 'Mtv_NMatch14';  % 14 faces (except for 1 and 16) while matching the number of presentation of each face   (See Model_Mtv_NMatch14.m)
%design_name = 'Mtv_AllFaces14';  % all presentations of 14 faces (except for 1 and 16) including the F0,F1, anf F2 events       (See Model_Mtv_AllFaces14.m)
% 2. contrasts names (which faces will be included to make an RDM matrix)
cont_names={'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8',...
                        'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15'};
%. 3. select the process to perform
modelspec=1; %0 if you want to generate new contrasts based on what has been modeled before, otherwise 1
contrstgen=1; %0 if you don't need to make con files from the model otherwise 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic model specification
for s = 1:nsubj
    if modelspec
        disp(['%%  Starting model ',design_name,' for subject : ', subj{s},'   %%']);
        if ~exist([an_path, design_name, fs, subj{s}, fs],'dir')
            mkdir([an_path, design_name, fs, subj{s}, fs]) 
        end
        
        fmri_spec.dir             = {[an_path,design_name,fs,subj{s},fs]}; % directory for SPM.mat file
        fmri_spec.timing.units    = 'secs';
        fmri_spec.timing.RT       = 1.3;        % TR
        fmri_spec.timing.fmri_t   = 46;        % microtime resolution (time-bins per scan) = slices
        fmri_spec.timing.fmri_t0  = 23;       % reference slice from pre-processing
        fmri_spec.mask            = {'/usr/local/MATLAB/spm12/tpm/mask_ICV.nii'};
        fmri_spec.mthresh         = -Inf;        %Mask
        fmri_spec.volt            = 1;

    % number of scans and data
        for sess = 1:nblocks
            epiimgs =spm_select('List', [data_path,subj{s}, fs, info.prefix.Run,num2str(sess)], '^wua.*.img');  %collecting epi images before smoothing (not swua.*.img)
            epiimgs=strcat([data_path,subj{s}, fs, info.prefix.Run,num2str(sess),fs], epiimgs);
            for epi=1:length(epiimgs)
                fmri_spec.sess(sess).scans{epi, 1}=epiimgs(epi, :);     
            end
        end

    %% Loading mat file extracted from behavior results (it has the matrix called session(sess))
        load ([ProjSet.Metapath, fs, 'NoPhi', fs,'Onset_phi_', subj{s}, '.mat']);
        session=Session.se;
        
    %% Calling Model specification function
       eval(['[fmri_spec, conditions] = Model_' design_name '(subj{s}, nblocks, session, fmri_spec);']);

    %% Run design specification
        matlabbatch{1}.spm.stats.fmri_spec = fmri_spec;
        disp(['%%%%%%%%% Starting design spec for subj: ', subj{s},'%%%%%%%%%']);
        batchname=strcat(fmri_spec.dir{1},design_name,'_model.mat');
        save(batchname,'matlabbatch');
        modelconf{s,1}=matlabbatch;
        clear matlabbatch;
        clear fmri_spec;
      end
end

%% Turn on the parallel computing
    pflag = gcp('nocreate');
    if isempty(pflag)
        poolsize=0;
    else
        poolsize=pflag.NumWorkers;
    end
    
%% Estimate
if modelspec
    Respath=fullfile(an_path, design_name); % where Beta files will be saved
    save([Respath, fs, 'MDconfg_', datestr(now),'.mat'], 'modelconf'); %Save batch file
    parfor si=1:length(subj)
        matlabbatch=modelconf{si,1};
        spm_jobman('run',matlabbatch);
        disp(subj{si});
    end
    clear matlabbatch;
    for sj=1:length(subj)
        estcong{sj,1}.spm.stats.fmri_est.spmmat = {[Respath, fs, subj{sj}, fs, 'SPM.mat']};
        estcong{sj,1}.spm.stats.fmri_est.write_residuals = 0;
        estcong{sj,1}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{sj}=estcong{sj,1};
    end
    spm_jobman('run',matlabbatch);
end

%% Contrast generation  
if contrstgen
     for s = 1:nsubj
%% Generate separate con and spmT files per block (session).
    % Make a folder for the contrast results if not made
        if ~exist([an_path, design_name, fs, subj{s}, fs],'dir')
            mkdir([an_path, design_name, fs, subj{s}, fs])
        end
        spmpath=[an_path, design_name, fs, subj{s}, fs];

        %load in lvl1 mat file
        spm.stats.con.spmmat = {fullfile(spmpath, 'SPM.mat')};
        load( fullfile(spmpath, [design_name '_model.mat']) );
        fmri_spec=matlabbatch{1,1}.spm.stats.fmri_spec;
        
        % Sanity check : the contrasts should be defined identically across
        % blocks (sessions), which is important for computing dissimialrity
        % across blocks (sessions)
        for bl=1:nblocks
            contname{bl}={fmri_spec.sess(bl).cond.name};
            if bl~=nblocks && nblocks~=1
                if mean(strcmp({fmri_spec.sess(bl).cond.name}, {fmri_spec.sess(bl+1).cond.name}))~=1
                    flag=1;
                end
            elseif bl==nblocks && nblocks~=1
                if mean(strcmp({fmri_spec.sess(1).cond.name}, {fmri_spec.sess(2).cond.name}))~=1
                    flag=1;
                end
            end
        end
        if flag
            error('Contrasts are different across blocks');
        else
            cont_names = contname{1};
        end
        
        % Assign contrasts
        for ci=1:numel(cont_names)
            extracont=numel(matlabbatch{1,1}.spm.stats.fmri_spec.sess(1).cond)-numel(cont_names); % in case there are regressors of non-interests are modeled together.
            spm.stats.con.consess{ci}.tcon.name =cont_names{ci};
            cntw=zeros(1, numel(cont_names)+motion_reg+extracont);
            cntw(1,ci)=1;
            spm.stats.con.consess{ci}.tcon.weights =cntw;
            spm.stats.con.consess{ci}.tcon.sessrep = 'sess'; % Generate separate con and spmT files per block(session).  put 'repl' to generate a combined con and spmT across blocks
        end
        spm.stats.con.delete = 1;

	%% Run Contrasts generate Batch
        disp(['%%%%%%%%% Running Contrasts for subj: ', subj{s},'%%%%%%%%%']);

        matlabbatch{1}.spm=spm;
        batchname=strcat(fmri_spec.dir{1},design_name,'_contrast.mat');
        save(batchname,'matlabbatch');
        spm_jobman('run',matlabbatch);
     end
     
     %% save GLM specification
        savepath=[an_path, design_name, fs];
        glmname=['GLM_' design_name '.mat'];
        if exist(glmname, 'file') ~= 2
            glmspec.npm_reg=cont_names; %non-parametric regressors
            glmspec.cont=cont_names; %Regressor of interests (number is coresponding to the con files)
        save(fullfile(savepath, glmname), 'glmspec');
        end
    clear fmri_spec epiimgs session Session matlabbatch spm SPM batchname;
end