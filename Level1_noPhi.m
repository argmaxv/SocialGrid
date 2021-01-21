%% 1st level analysis with generating contrasts
% GLM1 F1 and F2 activity is modulated by cos 6 theta
% by A PARK 12/28/2017

%% Path setup
clear;
close all;
clc;
    
[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[subj0, nsubj0] = CallSubj_PS;

nblocks = info.Nday*info.Nses;
base_path=ProjSet.basepath;     
data_path =ProjSet.DATApath;   
an_path =[ProjSet.ANA1path, 'NoPhi', fs]; %Where results will be saved (the Beta files)
Behavior_path=ProjSet.Metapath;             %  Where onset files are saved (from text2mat function)
model_path=ProjSet.MODELpath;              % Where model and contrast specification functions are
addpath (model_path);
periodicity=info.periodicity;                         % 6 folds periodicity

%%%%%%%%%%%%%%%%%%% Input the following information %%%%%%%%%%%%%%%%%%%%%%%%%
design_name = ['Grid', num2str(periodicity), '_noPhi']; %Save folder name
model_name = 'NoPhi_5s'; % Call Model_NoPhi_5s.m for model specification and Contrast_NoPhi_5s.m for contrast specification
cont_names={'CvF01F02', 'SvF01F02', 'ftest_CvSvF01F02'}; % Regressors of interests
modelspec=0; %0 if you want to generate new contrasts based on what has been modeled before, otherwise 1
contrstgen=1; %0 if you don't need to make con files from the model otherwise 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic model specification
for s = 1:nsubj0

    if modelspec
        disp(['%%  Starting model ',design_name,' for subject : ', subj0{s},'   %%']);
        if ~exist([an_path, design_name, fs, subj0{s}, fs],'dir')
            mkdir([an_path, design_name, fs, subj0{s}, fs]) 
        end
        
        fmri_spec.dir             = {[an_path,design_name,fs,subj0{s},fs]}; % directory for SPM.mat file
        fmri_spec.timing.units    = 'secs';
        fmri_spec.timing.RT       = 1.3;        % TR
        fmri_spec.timing.fmri_t   = 46;        % microtime resolution (time-bins per scan) = slices
        fmri_spec.timing.fmri_t0  = 23;       % reference slice from pre-processing
        fmri_spec.mask            = {[ProjSet.spmdir, fs, 'tpm', fs, 'mask_ICV.nii']};
        fmri_spec.mthresh         = -Inf; %MASK .2; %.8 or -Inf
        fmri_spec.volt            = 1;

       %% Get raw image data from the Data folder
            % integrating Day1 and Day2 scans into a single GLM
            
        datapath{1}=[data_path,[info.prefix.day1, subj0{s}]];   %Day1 scans
        datapath{2}=[data_path,[info.prefix.day2, subj0{s}]];   %Day2 scans
        for rp=1:2
            if ~exist(datapath{rp},'dir')
                error([datapath{rp} ' is not exist']);
            end
        end

        for sess = 1:nblocks % while sess varys in a range of 1 to 6, newses indicate what is the current session in each day scan (1~3)
            if sess<(nblocks/2)+1
                subj{s}={[info.prefix.day1, subj0{s}]};
                newses=sess;
            else
                subj{s}={[info.prefix.day2, subj0{s}]};
                newses=sess - (nblocks/2);
            end
            epiimgs =spm_select('List', strcat(data_path, subj{s}, fs, info.prefix.Run, num2str(newses)), '^swua.*.img');
            epiimgs=strcat(data_path, subj{s}, fs, info.prefix.Run, num2str(newses), fs, epiimgs);
            for epi=1:length(epiimgs)
                fmri_spec.sess(sess).scans{epi, 1}=epiimgs{epi, :};     
            end
        end

    %% Get the behavioral data (from each day's mat file)
        load ([Behavior_path, 'NoPhi', fs, 'Onset_phi_', info.prefix.day1, subj0{s}, '.mat']);
        sessionDay1 = Session.se;
        load ([Behavior_path, 'NoPhi', fs, 'Onset_phi_', info.prefix.day2, subj0{s}, '.mat']);
        sessionDay2 = Session.se;
    %% Calling Model specification function
        % Details of the GLM specification are determined in Model_[model_name].m
       eval(['[fmri_spec, conditions] = Model_' model_name '(subj0{s}, nblocks, sessionDay1, sessionDay2, fmri_spec, periodicity);']);

   %% Run design specification
        matlabbatch{1}.spm.stats.fmri_spec = fmri_spec;
        disp(['%%%%%%%%% Starting design spec for subj: ', subj0{s},'%%%%%%%%%']);
        disp(['Model spec: periodicity :', num2str(periodicity), ' Phi = 0 rad']);
        batchname=strcat(fmri_spec.dir{1},design_name,'_model.mat');        % Batch file for the GLM specification is saved
        save(batchname,'matlabbatch');                
        spm_jobman('run',matlabbatch);
        clear matlabbatch; clear fmri_spec; 

    %% Estimate
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[an_path,design_name,'/',subj0{s},'/SPM.mat']};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        disp(['%%%%%%%%% Starting estimation for subj: ', subj0{s},' %%%%%%%%%']);
        spm_jobman('run',matlabbatch);

    end

%% Contrast generation
        if contrstgen
            %Make a folder for the contrast results
            if ~exist([an_path, design_name, fs, subj0{s}, fs],'dir')
                mkdir([an_path, design_name, fs, subj0{s}, fs])
            end

            %load in lvl1 mat file
            spm.stats.con.spmmat = {[an_path design_name fs subj0{s} fs 'SPM.mat']};
            load([an_path, design_name, fs, subj0{s}, fs, design_name,'_model.mat']);

            % make the array of condition names    
            cntm={}; cntmx={}; im=1; imx=1;
            fmri_spec=matlabbatch{1,1}.spm.stats.fmri_spec;

            %replicate contrasts to other blocks
            for cntidx=1:length(fmri_spec.sess(1).cond) 
                cntm{im}=fmri_spec.sess(1).cond(cntidx).name;
                im=im+1; %non-parametric regressors
                cntmx{imx}=fmri_spec.sess(1).cond(cntidx).name;
                imx=imx+1;
                if ~isempty(fmri_spec.sess(1).cond(cntidx).pmod) %length(fmri_spec.sess(1).cond(cntidx).pmod)>0
                    for cntidy=1:length(fmri_spec.sess(1).cond(cntidx).pmod)
                        cntmx{imx}=fmri_spec.sess(1).cond(cntidx).pmod(cntidy).name; %non + parametric regressors
                        imx=imx+1; 
                    end
                else
                end
            end

        % Calling contrast generate function
        eval(['[spm] = Cont_' model_name '(cont_names, cntmx, spm);']); % contrasts are determined by the regressors of interests

        % Run Contrasts generate Batch
        disp(['%%%%%%%%% Running Contrasts for subj: ', subj0{s},'%%%%%%%%%']);
        disp(['Contrasts generations: periodicity :', num2str(periodicity), ' Phi = 0 rad']);

        matlabbatch{1}.spm=spm;
        spm_jobman('run',matlabbatch);
        resPath=cell2mat(matlabbatch{1,1}.spm.stats.con.spmmat);
        resPath=resPath(1:end-7);
        clear matlabbatch spm Session session epi epiimags fmri_spec

    end % contrastsgen
end%subject loop

%% Save GLM specification (for 2nd level analysis and Beta extraction for Get_Phi.m)
glmname=fullfile([an_path, design_name, fs], ['GLM_' design_name '.mat']);
if exist(glmname, 'file') ~= 2 && contrstgen
    glmspec.npm_reg=cntm; %non-parametric regressors
    glmspec.pm_reg=cntmx; %Regressors (non+parametric regressors)
    glmspec.cont=cont_names; %Regressor of interests (number is coresponding to the con files)
save(glmname, 'glmspec');
end
