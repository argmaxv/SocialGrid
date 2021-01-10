% Get_Phi.m
% Extract betas from the ROI and compute Phi of each of the participants
% by A Park 2018/03/20

%% Path setup
clear
close all;
clc;

[ProjSet, fs, info, ROI, fname]=Call_default_PS;
base_path=ProjSet.basepath; 
data_path =ProjSet.DATApath;
an_path =ProjSet.ANA1path;
if ~exist('circ_mean.m','file') % include the circular stat toolbox into your path
    addpath(ProjSet.circstat);
end
[subj, subn] = CallSubj_PS;
nblocks = info.Nses; % Number of blocks
nDay=info.Nday; % Number of days

% Define ROI
roitype = ROI.Grid{1};
mask=[ProjSet.ROIpath, roitype, fs, [roitype '_roi.mat']];

%% Main
for periodicity=info.periodicity:info.periodicity  %6-folds periodicity
    
    disp(['Periodicity:' num2str(periodicity)]);
    Phi0path = ['Grid', num2str(periodicity), '_noPhi'];
    design_name = ['Grid', num2str(periodicity), '_F01F02_5s'];
    dscript='NoPhi';                    % Model discription defined in the GLM (Lv1 analysis)
    cont_names={'vF01F02'};   % regressor of interest (Defined in contrast specification in GLM (Lv1.)

        for s=1:length(subj)
            fprintf('\n\n\n')
            disp(subj{s})
            load([an_path, dscript, fs, Phi0path, fs, Phi0path, '.mat']); %load GLM specification defined in 1st Lv.
            
%% Beta file extraction
    % read SPM.mat
    % get information about Beta files
    % make a maks which inform which beta is associalted
    % with the sin theta and cos theta (assigning a code such as 11,21,1,22) 
                    clear SPM
                    clear BetaInfo BetaInfoidx
                    load([an_path, dscript, fs, 'Grid6_noPhi', fs, subj{s}, fs, 'SPM.mat'])

%% Step1. Based on the how the GLM was generated. 
    % Assign a specific code to the Betas of interests 
    % Otherwise, we assigned the code '0'
                  
                    for nReg=1:size(SPM.Vbeta,2) % number of Betas files
                        BetaInfo{nReg}=SPM.Vbeta(nReg).descrip(29:end-6); % Get the event type from Beta infomation: F1 or F2 e.g.  'spm_spm:beta (0003) - Sn(1) F1F2xCvF01F02^1*bf(1)'; 29:end-6 indicate the event such as 'F1F2xCvF01F02^1'
                        
                            if strcmp(BetaInfo{nReg},  'F1F2xSvF01F02^1')  % F1 and F2 activity modeled by Sin(periodicity*theta)
                                BetaInfoidx(nReg)=1;
                            elseif strcmp(BetaInfo{nReg},  'F1F2xCvF01F02^1')  % F2 activity modeled by Cos(periodicity*theta)
                                BetaInfoidx(nReg)=2;
                            else
                                BetaInfoidx(nReg)=0;    %otherwise 0
                            end
                            
                    end                        

%% Step2. Because the Beta of sin and cos thetas are generated in each or
    % nblocks*nDay blocks (by the design of the GLM), assign the codes to the
    % other blocks as well.

                    for block=1:nblocks*nDay 
                        nBperBlock=(size(SPM.Vbeta,2)-(nblocks*nDay))/(nblocks*nDay); % number of Beta in each block
                        for tr=((block-1)*nBperBlock)+1:(block*nBperBlock) % Beta number (beta_00xx) in each block 
                            if BetaInfoidx(tr)~=0 % if not Betas of non-interests
                                clear data meanB rawB info
                                data=[an_path, 'Ground0', fs, Phi0path, fs, subj{s}, fs, 'beta_' sprintf('%04d',tr) '.nii']; % Betas of interests
                                [meanB, rawB, info] = extract_voxel_values(mask,data); % Extract the Beta values in the mask (ROI)            

                                if BetaInfoidx(tr)==1                                   %'F1F2xSvF01F02^1'
                                    sbj(s).vF01F02.ses(block).Sm=meanB;                         
                                    sbj(s).vF01F02.ses(block).SB=rawB.I.Ya;                             
                                elseif BetaInfoidx(tr)==2                             %'F1F2xCvF01F02^1'
                                    sbj(s).vF01F02.ses(block).Cm=meanB;                         
                                    sbj(s).vF01F02.ses(block).CB=rawB.I.Ya;                             
                                end

                            end %if non interests
                        end %for Beta
                    end % for block

                    for bcon=1:numel(cont_names)
                        voxblock=[]; % the variable which will contain in which block the beta was extracted.
                        voxphi=[]; % the variables which will contain the phi of each voxels in the ROI.

%% Step3. Compute the phi [atan(SinBeta/CosBeta)]/periodicity
                        for block=1:nblocks*nDay    
                            sbj(s).(cont_names{bcon}).ses(block).Phi=atan(sbj(s).(cont_names{bcon}).ses(block).SB./sbj(s).(cont_names{bcon}).ses(block).CB);                                                                               % 1.2 and 2.2 (prefered) - 3.2                            
                            phi_sample=sbj(s).(cont_names{bcon}).ses(block).Phi(~isnan(sbj(s).(cont_names{bcon}).ses(block).Phi)); 
                            sbj(s).(cont_names{bcon}).ses(block).phi_rmmissing=phi_sample;
                            voxblock=[voxblock; block*ones(1,length(phi_sample))];
                            voxphi=[voxphi; phi_sample];

        % 3.0 All mean (not cross-validation)
                            [mu, ul, ll]=circ_mean(phi_sample');
                            sbj(s).(cont_names{bcon}).ses(block).Phi_CirMean=mu/periodicity;
                            sbj(s).(cont_names{bcon}).ses(block).Phi_Cir95CI=[ul, ll];
                            sbj(s).(cont_names{bcon}).ses(block).Phi_CirSTD=circ_std(phi_sample');
                            clear phi_sample mu ul ll muw ulw llw;
                        end
                        
        % 3.1 CVxB -> cross block in the same day

                    for nd=1:nDay
                        for block=1:nblocks
                            clear phi_sample phi_w curBlock
                            curBlock=(nd-1)*nblocks + block;
                            xBlockidx0=[1:nblocks]+(nd-1)*nblocks; % Blocks acquired in the same day scan: 1,2,3 or 4,5,6
                            xBlockidx=xBlockidx0(~ismember(xBlockidx0, curBlock)); % Creating Cross-validation xBlocks in the same day e.g. when current block is 1 -> return 2,3 ... when 6 -> 4,5
                            phi_sample=voxphi(ismember(voxblock, xBlockidx)); % select the phi samples from differnt blocks acquired in the same day

                            [mu, ul, ll]=circ_mean(phi_sample);
                            sbj(s).(cont_names{bcon}).ses(curBlock).Phi_xBCirMean=mu/periodicity;
                            sbj(s).(cont_names{bcon}).ses(curBlock).Phi_xBCir95CI=[ul, ll];
                            sbj(s).(cont_names{bcon}).ses(curBlock).Phi_xBCirSTD=circ_std(phi_sample);
                            clear mu ul ll muw ulw llw;
                        end
                    end
                        
        % 3.2 CVxD -> cross different days scans

                    xDayidx=[nblocks+1:nblocks*nDay;1:nblocks];
                    for nd=1:nDay
                        clear voxphiOtherDay voxphiwOtherDay
                        voxphiOtherDay=voxphi(xDayidx(nd,:),:);
                        voxblockOtherDay=voxblock(xDayidx(nd,:),:);
                        clear voxphiPerDay voxphiwPerDay voxblockPerDay
                        voxphiPerDay=voxphi(1+((nd-1)*nblocks):nblocks+((nd-1)*nblocks),:);
                        voxblockPerDay=voxblock(1+((nd-1)*nblocks):nblocks+((nd-1)*nblocks),:);
                        for block=1:nblocks
                            clear samplephi_cv mucv ulcv llcv curBlock;
                            curBlock=(nd-1)*nblocks + block;
                            [muxd, ulxd, llxd]=circ_mean(voxphiOtherDay(:));
                            sbj(s).(cont_names{bcon}).ses(curBlock).Phi_xDCirMean = muxd/periodicity;
                            sbj(s).(cont_names{bcon}).ses(curBlock).Phi_xDCir95CI = [ulxd llxd];
                            sbj(s).(cont_names{bcon}).ses(curBlock).Phi_xDCirSTD=circ_std(voxphiOtherDay(:));
                        end
                    end                       
            end
    end %for
    eval([design_name '=sbj']);
    svpath = [ROIpath, fs, roitype];
    svfilename = ['Phi_', design_name, '_', roitype, '.mat'];
    if ~exist(svpath,'dir')
        mkdir(svpath);
    end
    save(fullfile(svpath, svfilename), design_name);
    clear sbj;
end % periodicity
