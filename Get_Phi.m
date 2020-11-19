% Get_Phi.m
% Extract betas from the ROI and compute Phi of each participants
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
nblocks = info.Nses;
nDay=info.Nday; % Number of days

% Define ROI
maskpath='/home/ravic/Pr4.PF/Imaging/Analysis/Analysis Phi/';
roitype = ['lERC_f100'];
ROIpath=[maskpath, 'GeneratedModel/'];
mask=[ROIpath, 'lERC_f100_roi.mat'];
dscript='NoPhi'; % model discription from the GLM
%%
for periodicity=info.periodicity:info.periodicity  %6-folds periodicity
    
    disp(['Periodicity:' num2str(periodicity)]);
    Phi0path = ['Grid', num2str(periodicity), '_noPhi'];
    design_name = ['Grid', num2str(periodicity), '_F01_F02_F12_FDec_5s'];
    %event_names={'vF01', 'vF02'}; % activity to model
    cont_names={'vF01F02'}; % regressor of interest

        for s=1:length(subj)
            fprintf('\n\n\n')
            disp(subj{s})
            %load([an_path, 'Ground0', fs, Phi0path, fs, 'GLM_Grid6_noPhi.mat']);
            load([an_path, dscript, fs, Phi0path, fs, 'GLM_Grid6_noPhi.mat']);
            
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
                        %if ~isempty(BetaInfo{nReg})
%                             if strcmp(BetaInfo{nReg}, 'F1') 
%                                 BetaInfoidx(nReg)=10;
%                             elseif strcmp(BetaInfo{nReg}, 'F2') 
%                                 BetaInfoidx(nReg)=20;

%                             if strcmp(BetaInfo{nReg}, 'F1xSvPF01^1')  % F1 activity modeled by Sin periodicity*theta F0F1 vector
%                                 BetaInfoidx(nReg)=11;
%                             elseif strcmp(BetaInfo{nReg}, 'F2xSvPF02^1')  % F2 activity modeled by Sin periodicity*theta F0F2 vector
%                                 BetaInfoidx(nReg)=21;
%                             elseif strcmp(BetaInfo{nReg}, 'F1xCvPF01^1')  % F1 activity modeled by Cos periodicity*theta F0F1 vector
%                                 BetaInfoidx(nReg)=12;
%                             elseif strcmp(BetaInfo{nReg}, 'F2xCvPF02^1')  % F2 activity modeled by Cos periodicity*theta F0F2 vector
%                                  BetaInfoidx(nReg)=22;
                        
                            if strcmp(BetaInfo{nReg},  'F1F2xSvF01F02^1')  % F1 and F2 activity modeled by Sin(periodicity*theta)
                                BetaInfoidx(nReg)=1;
                            elseif strcmp(BetaInfo{nReg},  'F1F2xCvF01F02^1')  % F2 activity modeled by Cos(periodicity*theta)
                                BetaInfoidx(nReg)=2;
                                
                                 
%                             elseif strcmp(BetaInfo{nReg}, 'F2xSvPFDec^1')  % F2 activity modeled by Sin periodicity*theta F0F2 vector
%                                 BetaInfoidx(nReg)=23;
%                             elseif strcmp(BetaInfo{nReg}, 'F2xCvPFDec^1')  % F2 activity modeled by Sin periodicity*theta F0F2 vector
%                                 BetaInfoidx(nReg)=24;
                            else %Btn
                                BetaInfoidx(nReg)=0;    %otherwise 0
                            end
                        %end
                    end
%                         for bcon=1:numel(cont_names)
%                             if bcon>(numel(cont_names)/2)
%                                 dayx=2;
%                             else
%                                 dayx=1;
%                             end
%                             %Cbetas(bcon)=find(strcmp(glmspec.pm_reg,[day_names{dayx}, 'C', cont_names{bcon}])); % number of beta_000x for Cos vec in the session 1
%                             %Sbetas(bcon)=find(strcmp(glmspec.pm_reg,[day_names{dayx}, 'S', cont_names{bcon}])); % number of beta_000x for Sin vec in the session 1
%                             Cbetas(bcon)=find(strcmp(glmspec.cont,[day_names{dayx}, 'C', cont_names{bcon}])); % number of beta_000x for Cos vec in the session 1
%                             Sbetas(bcon)=find(strcmp(glmspec.cont,[day_names{dayx}, 'S', cont_names{bcon}])); % number of beta_000x for Sin vec in the session 1
%                             
%                             find(strcmp(glmspec.pm_reg,'CvPF01'));
%                             
%                             glmspec.pm_reg
%                         end
                        
%                         clear Cbetas Sbetas
%                         Cmtx0=[0,1,0,0,1,0,1,0,0];
%                         Smtx0=[0,0,1,0,0,1,0,1,0];
%                         F10=2;
%                         F20=4;
%                         
% %                         Cmtx0=repmat([[0,1,0,0,1,0,1,0,0], zeros(1,nummotion)],1,nblocks*nDay);
%                         Cbetas=find(Cmtx0==1);
% %                         Smtx0=repmat([[0,0,1,0,0,1,0,1,0], zeros(1,nummotion)],1,nblocks*nDay);
%                         Sbetas=find(Smtx0==1);
%                        Cmtx(1,:)=find(Cmtx0==1);
%                         dayidx=[];
%                         for ni=1:nDay
%                             dayidx=[dayidx, ones(1,numel(Cmtx(1,:))/nDay)*ni];
%                         end
%                         Cmtx(2,:)=dayidx;

                        
                        

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
                                        
                                    if BetaInfoidx(tr)==1
                                        sbj(s).vF01F02.ses(block).Sm=meanB;                         
                                        sbj(s).vF01F02.ses(block).SB=rawB.I.Ya;                             
                                        %sbj(s).vF01F02.ses(block).S_nBeta(1)=tr;
                                    elseif BetaInfoidx(tr)==2
                                        sbj(s).vF01F02.ses(block).Cm=meanB;                         
                                        sbj(s).vF01F02.ses(block).CB=rawB.I.Ya;                             
                                        %sbj(s).vF01F02.ses(block).C_nBeta(1)=tr;
                                    end
                                    
                                end %if non interests
                            end %for Beta
                        end % for block

                        for bcon=1:numel(cont_names)
                            voxblock=[]; % the variable which will contain in which block the beta was extracted.
                            voxphi=[]; % the variables which will contain the phi of each voxels in the ROI.
%% Step3. Compute the phi [tan-1(SinBeta/CosBeta)]/periodicity
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
% %                     %% figure
% %                     if drawfig==1
% %                         posx=0.02*s+0.05;
% %                         posy=0.2-0.02*s;
% %                         cir(s)=figure('Units', 'normal', 'Position',[posx posy .65 .85]);%[500 400 600 350]);
% % 
% %                         for bcon=1:numel(cont_names)
% %                             subplot(4,numel(cont_names),bcon);
% %                             for block=1:nblocks
% %                                 %polarhistogram(sbj(s).cont(bcon).ses(block).phi_sample*periodicity, 'DisplayStyle', 'stairs'); %*periodicity
% %                                 y=polarhistogram(sbj(s).cont(bcon).ses(block).phi_sample/periodicity, 'DisplayStyle', 'stairs'); %-30 ~+30
% %                                 %rose(sbj(s).cont(bcon).ses(block).phi_sample);
% %                                 hold on;
% %                             end
% %                             hold off;
% %                             title(cont_names{bcon});
% %                             
% %                             subplot(4,numel(cont_names),numel(cont_names)+bcon);
% %                             for block=1:nblocks*nDay
% %                                 circ_plot(sbj(s).cont(bcon).ses(block).Phiw_CirMean);
% %                                 hold on;
% %                             end
% %                             hold off;
% %                             title([cont_names{bcon},' wSelf']);
% %                             
% %                             subplot(4,numel(cont_names),2*numel(cont_names)+bcon);
% %                             for block=1:nblocks*nDay
% %                                 circ_plot(sbj(s).cont(bcon).ses(block).Phiw_xBCirMean);
% %                                 hold on;
% %                             end
% %                             hold off;
% %                             title([cont_names{bcon},' wxBlock']);
% % 
% %                             subplot(4,numel(cont_names),3*numel(cont_names)+bcon);
% %                             for block=1:nblocks*nDay
% %                                 circ_plot(sbj(s).cont(bcon).ses(block).Phiw_xDCirMean);
% %                                 hold on;
% %                             end
% %                             hold off;
% %                             title([cont_names{bcon},' wxDay']);
% %                         end
% %                         
% %                         suptitle(subj{s});
% %                         
% %                         figpath=[ROIpath, 'Phi', fs, 'Period', num2str(periodicity), fs];
% % %                         if ~exist(figpath,'dir')
% % %                             mkdir(figpath);
% % %                         end
% %                         if ~exist([figpath, 'Phi', fs],'dir')
% %                             mkdir([figpath, 'Phi', fs]);
% %                         end
% %                         saveas(cir(s),[figpath, 'Phi', fs, subj{s}, '.eps']);
% %                         disp(['Figures are saved in ', figpath, 'Phi']);
% %                     end
% %                     
% % %             elseif strcmp(whattyp, 'spmF')
% % %                     %% spmF file extraction
% % %                     vectors={'vF01', 'vF02', 'vF01F02', 'vF12', 'vFDec'};
% % %                     %vectors={'vF01', 'vF02', 'vF01F02', 'vFDec'};
% % %                     for fcon=1:numel(vectors)
% % %                         disp(['       -----    ', vectors{fcon}]);
% % %                         data=[an_path, design_name, fs, subj{s}, fs, 'spmF_' sprintf('%04d',find(strcmp(glmspec.cont, ['C', vectors{fcon}]))) '.nii'];
% % %                         [Cm CR info] = extract_voxel_values(mask,data);
% % %                         data=[an_path, design_name, fs, subj{s}, fs, 'spmF_' sprintf('%04d',find(strcmp(glmspec.cont, ['S', vectors{fcon}]))) '.nii'];
% % %                         [Sm SR info] = extract_voxel_values(mask,data);
% % % 
% % %                         cur.Phi=atan(SR.I.Ya./CR.I.Ya)/periodicity;
% % %                         cur.meanPhi=mean(cur.Phi(~isnan(cur.Phi)));
% % %                         cur.sePhi=std(cur.Phi(~isnan(cur.Phi)))/sqrt(length(cur.Phi));
% % %                         sbj(s).(vectors{fcon})=cur;
% % %                         clear cur;
% % %                         sbj(s).(vectors{fcon}).mni=SR.I.mni;
% % %                         sbj(s).(vectors{fcon}).xyz=SR.I.xyz;    
% % %                     end
% % %             else
% % %                         msg = 'Error: %% define whattyp %%';
% % %                         error(msg)
% % 
% % %             end %if
        end %for
        
        
        eval([design_name '=sbj']);
        svpath = [ROIpath, fs, roitype];
        svfilename = ['Phi_', design_name, '_', roitype, '.mat'];
        if ~exist(svpath,'dir')
            mkdir(svpath);
        end
        %save([ROIpath, fs, roitype(2:end) fs whattyp, '_', design_name, roitype, nowis, '.mat'], design_name, 'subj');
        save(fullfile(svpath, svfilename), design_name, 'subj');
        clear sbj;
%    end %roix
end % periodicity
