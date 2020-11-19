% 1st level analysis with generating contrasts
%
% by SPARK 12/28/2017
% make sure 'design_name' and have the model file in the folder, 'Model'
% nblocks=number of session

clear all;
clc;

  for periodicity=6:6
%       if periodicity > 9 % if you want to exclude a specific periodicity 
%       else

        %% Input information

	subj = {'F01', 'F02', 'F03', 'F04', 'F05', ...
            'F06', 'F07', 'F08', 'F09', 'F10', ...
            'F11', 'F12', 'F13', 'F14', 'F15', ...
            'F16', 'F17', 'F18', 'F19', 'F20'};
    
        nblocks = 6;

        modelspec=0; %0 if you want to generate new contrasts based on what has been modeled before, otherwise 1
        contrstgen=1; %0 if you don't need to make con files from the model otherwise 1

        %% Contrasts (regressors of interests)

        design_name = ['Grid', num2str(periodicity), '_noPhi'];
        model_name = ['Grid0_F01_F02_FDec_5s'];
%         cont_names={'CvF01', 'CvF02', 'CvF01F02', 'CvFDec',...
%             'SvF01', 'SvF02', 'SvF01F02', 'SvFDec',...
%             'vF01', 'vF02', 'vFDec'};  
        cont_names={'vF01', 'vF02', 'vFDec'};  

        %% Path
        fs=filesep;
        move = 0;
        base_path='/mnt/datashare/Pr4.PF/Imaging/';
        data_path = '/mnt/datashare/Pr4.PF/Imaging/Data/';
        an_path =[base_path 'Analysis/Analysis Lv1/Ground0/'];
        ons_dir=[base_path 'Onset'];
        addpath '/usr/local/MATLAB/spm12/'
        addpath ([base_path fs 'Batch' fs]);
        addpath ([base_path fs 'Batch' fs 'Model' fs]);

        %% Basic model specification
        for s = 1:length(subj)
            if s==1 && periodicity==4
            else
            if modelspec==1
                        disp(['%%  Starting model ',design_name,' for subject : ', subj{s},'   %%']);
                if ~exist([an_path, design_name, fs, subj{s}, fs],'dir')
                    mkdir([an_path, design_name, fs, subj{s}, fs]) 
                end

                cd([an_path,design_name,fs,subj{s},fs]);

                fmri_spec.dir             = {[an_path,design_name,fs,subj{s},fs]}; % directory for SPM.mat file
                fmri_spec.timing.units    = 'secs';
                fmri_spec.timing.RT       = 1.3;        % TR
                fmri_spec.timing.fmri_t   = 46;        % microtime resolution (time-bins per scan) = slices
                fmri_spec.timing.fmri_t0  = 23;       % reference slice from pre-processing
                fmri_spec.mask            = {'/usr/local/MATLAB/spm12/tpm/mask_ICV.nii'};% YOU CAN EITHER GENERATE YOUR OWN MASK OR USE THE ONE IN SPM
                fmri_spec.mthresh         = -Inf; %MASK .2; %.8 or -Inf
                fmri_spec.volt            = 1;

            % number of scans and data
                for sess = 1:nblocks
                    epiimgs =spm_select('List', [data_path,subj{s}, fs, 'R',num2str(sess)], 'swua.*.img');  
                    epiimgs=strcat([data_path,subj{s}, fs, 'R',num2str(sess),'/'], epiimgs);
                    for epi=1:length(epiimgs)
                        fmri_spec.sess(sess).scans{epi, 1}=epiimgs(epi, :);     
                    end
                end

            %% Loading mat file extracted from behavior results (it has the matrix called session(sess))
                load ([ons_dir, fs, 'fMRI_Bhv_mat', fs, 'Ground0', fs, 'Onset_phi3_', subj{s}, '.mat']);
                session=Session.se;

            %% Calling Model specification function
               eval(['[fmri_spec, conditions] = Model_' model_name '(subj{s}, nblocks, session, fmri_spec, periodicity);']);

            %% Run design specification
                matlabbatch{1}.spm.stats.fmri_spec = fmri_spec;
                disp(['%%%%%%%%% Starting design spec for subj: ', subj{s},'%%%%%%%%%']);
                disp(['Model spec: periodicity :', num2str(periodicity), ' Phi = 0 rad']);

                spm_jobman('run',matlabbatch);
                batchname=strcat(fmri_spec.dir{1},design_name,'_model.mat');
                %save([fmri_spec.dir{1},design_name,'_model.mat'],'matlabbatch');
                %save(batchname,'matlabbatch');
                save(batchname,'fmri_spec');
                clear matlabbatch; clear fmri_spec; 

                %load ('SPM.mat')

                %% Estimate
                matlabbatch{1}.spm.stats.fmri_est.spmmat = {[an_path,design_name,'/',subj{s},'/SPM.mat']};
                matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

                disp(['%%%%%%%%% Starting estimation for subj: ', subj{s},' %%%%%%%%%']);

                spm_jobman('run',matlabbatch);

                %save([subj{s} 'struct.mat'], '-struct', 'spm')

                %clear epiimgs fmri_spec R mot_reg;
                %clear matlabbatch;
                %save (sprintf('%s%s/%s/cond22.mat',an_path,design_name,subj{s}),'cond')
            end

            if contrstgen==1
                %% Contrast generation

                %Make a folder for the contrast results if not made
                    if ~exist([an_path, design_name, fs, subj{s}, fs],'dir')
                        mkdir([an_path, design_name, fs, subj{s}, fs])
                    end

                %change to new directory for saving
                    cd([an_path, design_name, fs, subj{s}, fs])

                %load in lvl1 mat file
                    spm.stats.con.spmmat = {[an_path design_name fs subj{s} fs 'SPM.mat']};
                    load([design_name,'_model.mat']);
                    %load (batchname); %fmri_spec

                    % make the array of condition names    
                cntm={}; cntmx={}; im=1; imx=1;
                %fmri_spec=matlabbatch{1,1}.spm.stats.fmri_spec;
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
                eval(['[spm] = Cont_' model_name '(cont_names, cntmx, spm);']);

            %% Run Contrasts generate Batch
                disp(['%%%%%%%%% Running Contrasts for subj: ', subj{s},'%%%%%%%%%']);
                disp(['Contrasts generations: periodicity :', num2str(periodicity), ' Phi = 0 rad']);

                matlabbatch{1}.spm=spm;
                spm_jobman('run',matlabbatch);

                clear matlabbatch spm Session session epi epiimags fmri_spec
            end % contrastsgen
        end % special condition
        end%subject loop

        % save GLM specification
        cd([an_path, design_name, fs]);
        glmname=['GLM_' design_name '.mat'];
        if exist(glmname) ~= 2
            glmspec.npm_reg=cntm; %non-parametric regressors
            glmspec.pm_reg=cntmx; %Regressors (non+parametric regressors)
            glmspec.cont=cont_names; %Regressor of interests (number is coresponding to the con files)
        save(glmname, 'glmspec');
        end
      end    
