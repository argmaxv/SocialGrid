%% 1st level analysis with generating contrasts
% Univariate GLM2 (cos(6[theta-phi] only) + GLM3 (GP, Easiness only)
% Applying Phi estimated from EC roi with cross validation (across blocks)
% identidying the brain areas encoding cos(6[theta-phi], GP, Easiness, Euc
% by SPARK 12/28/2017

clear;
close all;
clc;

%% Path setup
[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[subj, nsubj] = CallSubj_PS;

%%%%%%%%%%%%%%%%%%% Input the following information %%%%%%%%%%%%%%%%%%%%%%%%%
ROIpath =ROI.Grid{1}; %EC
design_name = 'Grid_PhixB'; % Call Model and Contrast functions
% See ../Batch/Model/Model_Grid_PhixB.m for GLM design
% See ../Batch/Model/Contrast_Grid_PhixB.m for contrasts setup
cont_names={'vF0102', 'GP', 'Easiness', 'Euc'}; % Regressors of interests 
% cont_names={'vF0102'}; % GLM2
% cont_names={'GP', 'Easiness'}; % GLM3
modelspec=1; %0 if you want to generate new contrasts based on what has been modeled before, otherwise 1
contrstgen=1; %0 if you don't need to make con files from the model otherwise 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nblocks = info.Nday*info.Nses;
CrossValidationOption= 'PhixB'; % Phi across blocks
resfolder =[ROI.Grid{1}, '_', CrossValidationOption];
base_path=ProjSet.basepath;     
data_path =ProjSet.DATApath;   
an_path =[ProjSet.ANA1path];            % Where results of lv1 analysis are
Behavior_path=ProjSet.Metapath;      % Where onset files are saved (from text2mat function)
model_path=ProjSet.MODELpath;       % Where model and contrast specification functions are
addpath (model_path);
periodicity=info.periodicity;                  % 6 folds periodicity
Respath = [an_path, resfolder];           % Where results will be saved (the Beta files)

%% Basic model specification
if modelspec==1
        disp(['********* Applying Phi (cross blocks) estimated from the , ' ROIpath ' *********']);

%% Turn on the parallel computing
    pflag = gcp('nocreate');
    if isempty(pflag)
        poolsize=0;
    else
        poolsize=pflag.NumWorkers;
    end
    
	for s = 1:nsubj
        disp(['%%  Starting model ',design_name,' for subject : ', subj{s},'   %%']);
        if ~exist([Respath, fs, design_name, fs, subj{s}, fs],'dir')
            mkdir([Respath, fs, design_name, fs, subj{s}, fs]) 
        end
        cd([Respath, fs, design_name, fs, subj{s}, fs]);
        fmri_spec.dir             = {[Respath, fs, design_name, fs, subj{s}, fs]}; % directory for SPM.mat file
        fmri_spec.timing.units    = 'secs';
        fmri_spec.timing.RT       = 1.3;        % TR
        fmri_spec.timing.fmri_t   = 46;        % microtime resolution (time-bins per scan) = slices
        fmri_spec.timing.fmri_t0  = 23;       % reference slice from pre-processing
        fmri_spec.mask            = {'/usr/local/MATLAB/spm12/tpm/mask_ICV.nii'};
        fmri_spec.mthresh         = -Inf; 
        fmri_spec.volt            = 1;

        % number of scans and data
        for sess = 1:nblocks
            clear subjt;
            if sess<(nblocks/2)+1
                subjt{s}={[info.prefix.day1, subj{s}]};
                newses=sess;
            else
                subjt{s}={[info.prefix.day2, subj{s}]};
                newses=sess - (nblocks/2);
            end
            epiimgs =spm_select('List', strcat(data_path, subjt{s}, fs, info.prefix.Run, num2str(newses)), '^swua.*.img');  
            epiimgs=strcat(data_path, subjt{s}, fs, info.prefix.Run, num2str(newses), fs, epiimgs);
            for epi=1:length(epiimgs)
                fmri_spec.sess(sess).scans{epi, 1}=epiimgs{epi, :};     
            end
        end

    %% Loading mat file extracted from behavior results (it has the matrix called session(sess))
        onsettype = 'Onset_phi_';
        load ([Behavior_path, ROIpath, fs, CrossValidationOption, fs, onsettype, info.prefix.day1, subj{s}, '.mat']);
        sessionDay1=Session.se;
        load ([Behavior_path, ROIpath, fs, CrossValidationOption, fs, onsettype, info.prefix.day2, subj{s}, '.mat']);
        sessionDay2=Session.se;
                
    %% Calling Model specification function
       eval(['[fmri_spec, conditions] = Model_' design_name '(subj{s}, nblocks, sessionDay1, sessionDay2, fmri_spec, periodicity);']);

    %% Run design specification
        matlabbatch{1}.spm.stats.fmri_spec = fmri_spec;
        disp(['%%%%%%%%% Starting design spec for subj: ', subj{s},'%%%%%%%%%']);
        batchname=strcat(fmri_spec.dir{1},design_name,'_model.mat');
        save(batchname,'matlabbatch');
        modelconf{s,1}=matlabbatch;
        clear matlabbatch; clear fmri_spec;
	end
    save([fullfile(Respath, design_name), fs, 'MDconfg_', datestr(now),'.mat'], 'modelconf');
end
        
if modelspec==1    
    parfor si=1:length(subj)
        matlabbatch=modelconf{si,1};
        spm_jobman('run',matlabbatch);
    end

   %% Estimate
    clear matlabbatch;
    for sj=1:length(subj)
        estcong{sj,1}.spm.stats.fmri_est.spmmat = {[Respath, fs, design_name, fs, subj{sj}, fs, 'SPM.mat']};
        estcong{sj,1}.spm.stats.fmri_est.write_residuals = 0;
        estcong{sj,1}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{sj}=estcong{sj,1};
    end
    spm_jobman('run',matlabbatch);
end        

%% Contrast generation
if contrstgen==1
	for s = 1:length(subj)
        %Make a folder for the contrast results if not made
        if ~exist([Respath, fs, design_name, fs, subj{s}, fs],'dir')
            mkdir([Respath, fs, design_name, fs, subj{s}, fs])
        end
        
        %Load in lvl1 mat file
        spm.stats.con.spmmat = {fullfile([Respath, fs, design_name fs subj{s} fs], 'SPM.mat')};
        load(fullfile([Respath, fs, design_name fs subj{s} fs], [design_name,'_model.mat']));
        cntm={}; cntmx={}; im=1; imx=1;
        fmri_spec=matlabbatch{1,1}.spm.stats.fmri_spec;
        for cntidx=1:length(fmri_spec.sess(1).cond) %if contrast of session 2 is different from that of session 2 then reprogram the following for loop.
            cntm{im}=fmri_spec.sess(1).cond(cntidx).name;
            im=im+1; %non-parametric regressors
            cntmx{imx}=fmri_spec.sess(1).cond(cntidx).name;
            imx=imx+1;
            if ~isempty(fmri_spec.sess(1).cond(cntidx).pmod) 
                for cntidy=1:length(fmri_spec.sess(1).cond(cntidx).pmod)
                    cntmx{imx}=fmri_spec.sess(1).cond(cntidx).pmod(cntidy).name; %non + parametric regressors
                    imx=imx+1; 
                end
            else
            end
        end

	% Calling contrast generate function
        eval(['[spm] = Cont_' design_name '(cont_names, cntmx, spm);']);

	%% Run Contrasts generate Batch
        disp(['%%%%%%%%% Running Contrasts for subj: ', subj{s},'%%%%%%%%%']);
        contconf{s,1}.spm=spm;
        clear matlabbatch spm Session session epi epiimags fmri_spec
        end %subject loop
    save([fullfile(an_path, resfolder, design_name), fs, cont_name, '_Cont.mat'], 'contconf');
    for sk=1:length(subj)
        clear matlabbatch
        matlabbatch{1}=contconf{sk,1};
        spm_jobman('run',matlabbatch);
    end %subject loop
end % contrastsgen

%% Save GLM specification
cd([Respath, fs, design_name, fs]);
glmname=['GLM_' design_name '.mat'];
if exist(glmname) ~= 2 && contrstgen
    glmspec.npm_reg=cntm; %non-parametric regressors
    glmspec.pm_reg=cntmx; %Regressors (non+parametric regressors)
    glmspec.cont=cont_names; %Regressor of interests (number is coresponding to the con files)
    save(glmname, 'glmspec');
end