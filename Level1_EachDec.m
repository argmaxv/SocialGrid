% 1st level analysis with generating contrasts
%
% by SPARK 12/28/2017
% make sure 'design_name' and have the model file in the folder, 'Model'
% attention! I select the ^wus images not the smoothed images.
% In the process of Contrast, it won't read any con description but, make T for all regressors except with those with motions/sessions. 
% nblocks=number of session

clear all;
clc;

%% Path setup
[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[subj, nsubj] = CallSubj_PS;

%%%%%%%%%%%%%%%%%%% Input the following information %%%%%%%%%%%%%%%%%%%%%%%%%
ROIpath =ROI.Grid{1}; %EC
design_name = {'EachDec'};%{'EachDec_GP'};

modelspec=1; %0 if you want to generate new contrasts based on what has been modeled before, otherwise 1
contrstgen=0; %0 if you don't need to make con files from the model otherwise 1
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
    
%% Turn on the parallel computing
    pflag = gcp('nocreate');
    if isempty(pflag)
        poolsize=0;
    else
        poolsize=pflag.NumWorkers;
    end
	dsgn=1;
%for dsgn=1:numel(design_name)
    %% Model specification
if modelspec ==1
        for s = 1:length(subj0)
            disp(['%%  Starting model ',design_name{dsgn},' for subject : ', subj0{s},'   %%']);
            if ~exist([an_path, design_name{dsgn}, fs, subj0{s}, fs],'dir')
                mkdir([an_path, design_name{dsgn}, fs, subj0{s}, fs]) 
            end

            fmri_spec.dir             = {[an_path,design_name{dsgn},fs,subj0{s},fs]}; % directory for SPM.mat file
            fmri_spec.timing.units    = 'secs';
            fmri_spec.timing.RT       = 1.3;        % TR
            fmri_spec.timing.fmri_t   = 46;        % microtime resolution (time-bins per scan) = slices
            fmri_spec.timing.fmri_t0  = 23;       % reference slice from pre-processing
            fmri_spec.mask            = {[ProjSet.spmdir, fs, 'tpm', fs, 'mask_ICV.nii']};
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
        
        % Loading mat file extracted from behavior results (it has the matrix called session(sess))
            onsettype = 'Onset_phi_';
            load ([Behavior_path, ROIpath, fs, CrossValidationOption, fs, onsettype, info.prefix.day1, subj{s}, '.mat']);
            sessionDay1=Session.se;
            load ([Behavior_path, ROIpath, fs, CrossValidationOption, fs, onsettype, info.prefix.day2, subj{s}, '.mat']);
            sessionDay2=Session.se;
        
        % Calling Model specification function
            eval(['[fmri_spec, conditions] = Model_' design_name '(subj{s}, nblocks, sessionDay1, sessionDay2, fmri_spec);']);

        % Run design specification

            matlabbatch{1}.spm.stats.fmri_spec = fmri_spec;
            disp(['%%%%%%%%% Starting design spec for subj: ', subj0{s},'%%%%%%%%%']);

            batchname=strcat(fmri_spec.dir{1},design_name{dsgn},'_model.mat');
            save(batchname,'matlabbatch');
            modelconf{s,1}=matlabbatch;
            clear matlabbatch; clear fmri_spec;

    end % for loop subject
batchname=fullfile(an_path, design_name, 'Model_AllSubj.mat');
save(batchname{:},'modelconf');
end % if modelspec == 1

%% Model specification
if modelspec==1
    fprintf('Model Specification\n');
    parfor si=1:length(subj0)
        spmpath=fullfile(modelconf{si, 1}{1, 1}.spm.stats.fmri_spec.dir,'SPM.mat');
        if ~exist(spmpath{:},'file')
            fprintf('%s\t', subj0{si});
            matlabbatch=modelconf{si,1};
            spm_jobman('run',matlabbatch);
        else
            fprintf('%s\t', subj0{si});
        end
    end
end
%% Estimation
if modelspec==1
    clear matlabbatch;
    fprintf('Estimate\n');
    for sj=1:length(subj0)
        fprintf('%s\t', subj0{sj});
        estcong{sj,1}.spm.stats.fmri_est.spmmat = fullfile(an_path, design_name, subj0{sj}, 'SPM.mat');
        estcong{sj,1}.spm.stats.fmri_est.write_residuals = 0;
        estcong{sj,1}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{sj}=estcong{sj,1};
    end
    spm_jobman('run',matlabbatch);
end

%% Contrast generation
if contrstgen==1
    for sj=1:length(subj0)
    %Make a folder for the contrast results if not made
        if ~exist([an_path, design_name{dsgn}, fs, subj0{sj}, fs],'dir')
            mkdir([an_path, design_name{dsgn}, fs, subj0{sj}, fs])
        end

    %Take the specified GLM
        spmpath=[an_path, design_name{dsgn}, fs, subj0{sj}, fs];
        spm.stats.con.spmmat = {[spmpath 'SPM.mat']};
        load([spmpath design_name{dsgn} '_model.mat']);
        fmri_spec=matlabbatch{1,1}.spm.stats.fmri_spec;
        motion_reg=6;
        clear cont_names;
        cont_names={fmri_spec.sess(1).cond.name};
        newci=0;
        for ses=1:nblocks
            for ci=1:numel(cont_names) %the last one is btn
                if strcmp(cont_names{ci}, 'Btn')
                else
                    newci=newci+1;
                    spm.stats.con.consess{newci}.tcon.name =['S', num2str(ses), cont_names{ci}];
                    cntw=zeros(1, nblocks*(numel(cont_names)+motion_reg));
                    cntw(1, (ses-1)*(numel(cont_names)+motion_reg)+ci )=1;                
                    spm.stats.con.consess{newci}.tcon.weights =cntw;
                    spm.stats.con.consess{newci}.tcon.sessrep = 'none';    %'repl'; %replicates over sessions
                    spm.stats.con.delete = 1;   % 1:Delete exsisting con files
                end
            end
        end
    %% Run Contrasts generate Batch
        disp(['%%%%%%%%% Running Contrasts for subj: ', subj0{sj},'%%%%%%%%%']);
        matlabbatch{1}.spm=spm;
        spm_jobman('run',matlabbatch);
    end % subj
end %contrstgen
