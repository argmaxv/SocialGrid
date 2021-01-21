% 1st level analysis with generating contrasts
%
% by SPARK 12/28/2017
% make sure 'design_name' and have the model file in the folder, 'Model'
% attention! I select the ^wus images not the smoothed images.
% In the process of Contrast, it won't read any con description but, make T for all regressors except with those with motions/sessions. 
% nblocks=number of session

clear; close all; clc;

%% Path setup
[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[subj0, nsubj] = CallSubj_PS;

%%%%%%%%%%%%%%%%%%% Input the following information %%%%%%%%%%%%%%%%%%%%%%%%%
ROIpath =ROI.Grid{1}; %EC
design_name = {'GPEucF1F2_OnOff'};%GLM4 controlling for Euc;
model_name = {'GPEucF1F2_OnOff'}; %GLM4 controlling for Euc;    
cont_names={'F12on', 'F12off', 'F12onoff_diff', 'GPon', 'GPoff', 'GPonoff_diff'};    
modelspec=1; %0 if you want to generate new contrasts based on what has been modeled before, otherwise 1
contrstgen=1; %0 if you don't need to make con files from the model otherwise 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nblocks = info.Nday*info.Nses;
CrossValidationOption= 'PhixB'; % Phi across blocks
base_path=ProjSet.basepath;     
data_path =ProjSet.DATApath;   
an_path =[ProjSet.ANA1path];            % Where results of lv1 analysis are
Behavior_path=ProjSet.Metapath;      % Where onset files are saved (from text2mat function)

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
    load('Info_OnOff.mat'); % whether each event is on or off-grid
        for s = 1:length(subj0)
            clear OnOff_sub
            OnOff_sub=subjects(s);           
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
            numscan_st=0;
            for sess = 1:nblocks
                if sess<(nblocks/2)+1
                    subj{s}={[info.prefix.day1, subj0{s}]};
                    newses=sess;
                else
                    subj{s}={[info.prefix.day2, subj0{s}]};
                    newses=sess - (nblocks/2);
                end
                clear icnt numscan numscan_end epiimgs resimgs
                epiimgs =spm_select('List', strcat(data_path, subj{s}, fs, info.prefix.Run, num2str(newses)), '^swua.*.img');
                epiimgs=strcat(data_path, subj{s}, fs,  info.prefix.Run, num2str(newses), fs, epiimgs);
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
            eval(['[fmri_spec, conditions] = Model_' model_name{dsgn} '(subj0{s}, nblocks, sessionDay1, sessionDay2, fmri_spec, OnOff_sub);']);

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
    %% Contrast generation
    for sj=1:length(subj0)
        %Make a folder for the contrast results if not made
            if ~exist([an_path, design_name{dsgn}, fs, subj0{sj}, fs],'dir')
                mkdir([an_path, design_name{dsgn}, fs, subj0{sj}, fs])
            end

        %change to new directory for saving
            cd([an_path, design_name{dsgn}, fs, subj0{sj}, fs])

        %load in lvl1 mat file
            spm.stats.con.spmmat = {[an_path design_name{dsgn} fs subj0{sj} fs 'SPM.mat']};
            load([design_name{dsgn},'_model.mat']);
            %load (batchname); %fmri_spec
            spmpath=[an_path, design_name{dsgn}, fs, subj0{sj}, fs];
            load([spmpath design_name{dsgn} '_model.mat']);
            fmri_spec=matlabbatch{1,1}.spm.stats.fmri_spec;
            % make the array of condition names    
            cntm={}; cntmx={}; im=1; imx=1;
            for cntidx=1:length(fmri_spec.sess(1).cond) %if contrast of session 2 is different from that of session 2 then reprogram the following for loop.

                cntm{im}=fmri_spec.sess(1).cond(cntidx).name;
                im=im+1; %non-parametric regressors

                cntmx{imx}=fmri_spec.sess(1).cond(cntidx).name;
                imx=imx+1; 

                if length(fmri_spec.sess(1).cond(cntidx).pmod)>0
                    for cntidy=1:length(fmri_spec.sess(1).cond(cntidx).pmod)
                        cntmx{imx}=fmri_spec.sess(1).cond(cntidx).pmod(cntidy).name; %non + parametric regressors
                        imx=imx+1; 
                    end
                else
                end
            end

    %% Contrasts 

    % Calling contrast generate function
        eval(['[spm] = Cont_' model_name{dsgn} '(cont_names, cntmx, spm);']);

    %% Run Contrasts generate Batch
        disp(['%%%%%%%%% Running Contrasts for subj: ', subj0{sj},'%%%%%%%%%']);

        matlabbatch{1}.spm=spm;
        spm_jobman('run',matlabbatch);

        clear matlabbatch spm Session session epi epiimags fmri_spec
    end 
end% contrastsgen

    % save GLM specification
    cd([an_path, design_name{dsgn}, fs]);
    glmname=['GLM_' design_name{dsgn} '.mat'];
    if exist(glmname) ~= 2
        glmspec.npm_reg=cntm; %non-parametric regressors
        glmspec.pm_reg=cntmx; %Regressors (non+parametric regressors)
        glmspec.cont=cont_names; %Regressor of interests (number is coresponding to the con files)
        save(glmname, 'glmspec');
    end