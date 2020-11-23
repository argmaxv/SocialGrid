%% RSA ROI
% 3. Compute Tau
% by SPARK 1.Oct.2018

clear;
close all;
clc;

%% Setting

[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[subj, nsubj] = CallSubj_PS;
raw_path=ProjSet.DATApath;% raw data path
maskpath=ProjSet.ROIpath; %ROI path
Anal1path=ProjSet.ANA1path; %1st level Analysis path to find SPM.mat
Respath=ProjSet.Respath; %save path
RDMpath=ProjSet.Respath;% RSA results path 
ROIs=ROI.Ana; %ROI list - HCr,HCl, ECr,ECl,M1r, and M1l
nRun=info.Nses; %3
nDay=info.Nday; %2
Nperm=info.Nperm; % num permutation, 1000 


%%%%%%%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%%%%% 
models={'Mtv_NMatch14', 'Mtv_AllFaces14'}; % Names of 1st level analysis path (See 'design_name' in Multivariate_Level1.m)
%'Mtv_NMatch14' - 14 faces (except for 1 and 16) while matching the number of presentation of each face
%'Mtv_AllFaces14' - All presentations of 14 faces (except for 1 and 16) including the F0,F1, anf F2 events
rng=2:15; % 2~15 individuals among 16 individuals in social hierarchy
svoption=1; % 1 if you want to save the results of RDM as mat file. 0 otherwise.
histon=0; % 1 to draw histogram, 0 otherwise
maskoption=0; % In case you want to remove the cells in specific Euclidean distance in the analysis
svfile='Results_TauA'; %savefile name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% EucRDM_factorize generates the model RDM (2-D, 2 2-D and their permutations)
[E_rdm, D1_rdm, D2_rdm, E_rdm_perm, D1_rdm_perm, D2_rdm_perm]=EucRDM_factorize(rng, Nperm);

% In case you want to remove the cells in specific Euclidean distance in the analysis
if maskoption 
    themask=ones(size(E_rdm));
    themask(E_rdm==1)=nan;
    surfix='_rmEuc1';
else
    themask=ones(size(E_rdm));
    surfix='';
end
for m=1:numel(models)
    fprintf('\n %s   ', models{m});
    % Load RDM generated from RSA_2NoiseCov.m
    svpath=[RDMpath, models{m}, fs, 'RSA', fs]; % Save path
    if exist(fullfile([RDMpath, models{m}, fs], fname.RSAroi),'file')
        load(fullfile([RDMpath, models{m}, fs], fname.RSAroi)); % Where the RDM is; Load RDM (cross blocks)
    else
        error('The RDM is missing');
    end
        
    for r=1:numel(ROIs)
    fprintf('\n %s   ', ROIs{r});

        for s=1:nsubj
            fprintf('%d   ', s);
            clear RDM0cross
            RDMroisub0=squeeze(RDM.(ROIs{r}).RDM(s,:,:));
            RDMroisub=(RDMroisub0+RDMroisub0')/2;
            % For the next analysis, RSA_4Dissim_by_DistBin.m
            RDMcross.(ROIs{r}).norm(:,:,s)=normalize(RDMroisub,'range'); 
            RDMcross.(ROIs{r}).mean(:,:,s)=RDMroisub;

    %% Rank correation of Euclidean distances and two 1-D distances in each of two social hierarchy dimensions       
            E_rdm_mask0=E_rdm.*themask;
            RDMroisub_mask0=RDMroisub.*themask;
            E_rdm_mask=selectriu(E_rdm_mask0); % selectriu: make the RDM symmetry and choose the one side (right-top)
            RDMroisub_mask=selectriu(RDMroisub_mask0);
            tauEach.Ec0(s)=rsa.stat.rankCorr_Kendall_taua(RDMroisub_mask, E_rdm_mask); % with pairwise Euclidean distance model
            tauEach.D1c0(s)=rsa.stat.rankCorr_Kendall_taua(selectriu(RDMroisub), selectriu(D1_rdm) ); % with rank distandce in the Dimension 1
            tauEach.D2c0(s)=rsa.stat.rankCorr_Kendall_taua(selectriu(RDMroisub), selectriu(D2_rdm) ); % with rank distandce in the Dimension 2

   %% For comparion (E vs. D1; E vs. D2).
        % Input E and either D1 or D2 into the same GLM to compete to explain the variances of neural RDM.
            D1betas=glmfit( [selectriu(E_rdm), selectriu(D1_rdm)], selectriu(RDMroisub) );
            D2betas=glmfit( [selectriu(E_rdm), selectriu(D2_rdm)], selectriu(RDMroisub) );
            tauEach.E1beta(s)=D1betas(2);
            tauEach.D1beta(s)=D1betas(3);
            tauEach.E2beta(s)=D2betas(2);
            tauEach.D2beta(s)=D2betas(3);

    %% Computing permuted baseline
        %  Compute the Kendall's tau based on the permuted baseline

            for np=1:Nperm 
                clear E_rdmperm_mask0
                E_rdmperm_mask0=E_rdm_perm{np}.*themask;
                E_rdmperm_mask=selectriu(E_rdmperm_mask0);
               
                tauEach.permE(s,np)=rsa.stat.rankCorr_Kendall_taua(RDMroisub_mask,E_rdmperm_mask);
                tauEach.permD1(s,np)=rsa.stat.rankCorr_Kendall_taua(selectriu(RDMroisub), selectriu(D1_rdm_perm{np}) );
                tauEach.permD2(s,np)=rsa.stat.rankCorr_Kendall_taua(selectriu(RDMroisub), selectriu(D2_rdm_perm{np}) );    
            end
            
    % Correct according to the permuated baseline
            tauEach.E(s)=tauEach.Ec0(s)-mean(tauEach.permE(s,:));
            tauEach.D1(s)=tauEach.D1c0(s)-mean(tauEach.permD1(s,:));
            tauEach.D2(s)=tauEach.D2c0(s)-mean(tauEach.permD2(s,:));

        end %for subject s

    %% One-sample Wilcoxon signrank test 
             [p.E(r), ~]=rsa.stat.signrank_onesided(tauEach.E);
             [p.D1(r), ~]=rsa.stat.signrank_onesided(tauEach.D1);
             [p.D2(r), ~]=rsa.stat.signrank_onesided(tauEach.D2);
    % for Figures (Boxplot)         
             tau.(ROIs{r}).meanD1=mean(tauEach.D1);
             tau.(ROIs{r}).meanD2=mean(tauEach.D2);
             tau.(ROIs{r}).meanE=mean(tauEach.E);
             tau.(ROIs{r}).D1=tauEach.D1';
             tau.(ROIs{r}).D2=tauEach.D2';
             tau.(ROIs{r}).E=tauEach.E';
             tau.(ROIs{r}).seD1=SEM(tauEach.D1');
             tau.(ROIs{r}).seD2=SEM(tauEach.D2');
             tau.(ROIs{r}).seE=SEM(tauEach.E');
             tau.(ROIs{r}).E1beta=tauEach.E1beta';
             tau.(ROIs{r}).E2beta=tauEach.E2beta';
             tau.(ROIs{r}).D1beta=tauEach.D1beta';
             tau.(ROIs{r}).D2beta=tauEach.D2beta';
             tau.(ROIs{r}).tauRaw=tauEach;

      %% Comparion method 1 (E vs. D1; E vs. D2).
        % Paired signrank test
             [p.ED1(r), h.ED1(r), statsED1]=signrank(tauEach.E', tauEach.D1');
             [p.ED2(r), h.ED2(r), statsED2]=signrank(tauEach.E', tauEach.D2');
             [p.D1D2(r), h.D1D2(r), statsD1D2]=signrank(tauEach.D1', tauEach.D2');      
             z.ED1(r)=statsED1.zval;
             z.ED2(r)=statsED2.zval;
             z.D1D2(r)=statsD1D2.zval;

      %% Comparion method 2 (E vs. D1; E vs. D2).
        % Paired t test (Beta); Input E and either D1 or D2 into the same GL to compete to explain the variances of neural RDM.
        [~, p_glm1, ~, stat_glm1]=ttest(tauEach.E1beta, tauEach.D1beta);
        p.ED1_glm(r)=p_glm1;
        t.ED1_glm(r)=stat_glm1.tstat;

        [~, p_glm2, ~, stat_glm2]=ttest(tauEach.E2beta, tauEach.D2beta);
        p.ED2_glm(r)=p_glm2;
        t.ED2_glm(r)=stat_glm2.tstat;
    end %for ROIs r

    %% Save results to mat
    tau.All.p=p;
    tau.All.h=h;
    tau.All.z=z;
    tau.All.t=t;
    if svoption
        save(fullfile(svpath, [svfile, surfix, '.mat']), 'tau','Nperm');
            if ~maskoption
                % input data for RSA_4Dissim_by_DistBin.m
                save(fullfile([RDMpath, models{m}, fs], fname.Darwroi), 'RDMcross');
            end
    end
    tau_allModels.(models{m})=tau;
    
    %% Figures
    tau.All.E=[]; tau.All.E_mean=[]; tau.All.E_se=[]; 
    tau.All.D1=[]; tau.All.D1_mean=[]; tau.All.D1_se=[]; 
    tau.All.D2=[]; tau.All.D2_mean=[]; tau.All.D2_se=[]; 
    tau.All.E1_beta=[]; tau.All.E2_beta=[]; 
    tau.All.D1_beta=[]; tau.All.D2_beta=[]; 

    for r=1:numel(ROIs)
        tau.All.E=[tau.All.E, tau.(ROIs{r}).E];
        tau.All.D1=[tau.All.D1, tau.(ROIs{r}).D1];
        tau.All.D2=[tau.All.D2, tau.(ROIs{r}).D2];

        tau.All.E_mean=[tau.All.E_mean, tau.(ROIs{r}).meanE];
        tau.All.D1_mean=[tau.All.D1_mean, tau.(ROIs{r}).meanD1];
        tau.All.D2_mean=[tau.All.D2_mean, tau.(ROIs{r}).meanD2];

        tau.All.E_se=[tau.All.E_se, tau.(ROIs{r}).seE];
        tau.All.D1_se=[tau.All.D1_se, tau.(ROIs{r}).seD1];
        tau.All.D2_se=[tau.All.D2_se, tau.(ROIs{r}).seD2]; 

        tau.All.E1_beta=[tau.All.E1_beta,  tau.(ROIs{r}).E1beta];
        tau.All.E2_beta=[tau.All.E2_beta,  tau.(ROIs{r}).E2beta];
        tau.All.D1_beta=[tau.All.D1_beta,  tau.(ROIs{r}).D1beta];
        tau.All.D2_beta=[tau.All.D2_beta,  tau.(ROIs{r}).D2beta];
    end

    clear fig
    % Fig1. TauA for 2D Euclidean RDM
    i=1; pos = [-.2+(i-1*.3) .2 .15 .20];
    fig(i)=figure;
    set(gcf, 'Color', ones(1,3), 'Units', 'Normalized', 'Position', pos);
    boxplot(tau.All.E);
    hold on;
    errorbar([1:numel(ROIs)], tau.All.E_mean, tau.All.E_se, '.');
    title('Tau_A for Euclidean RDM');
    xticklabels(ROIs);
    hold off;
    if svoption
        saveas(fig(i), fullfile(svpath,[svfile, '_E', surfix]), 'epsc');
    end

    if ~maskoption
	% Fig2. TauA for two 1D RDM
        i=i+1;  pos = [-.2+(i-1*.3) .2 .35 .20];
        fig(i)=figure;
        set(gcf, 'Color', ones(1,3), 'Units', 'Normalized', 'Position', pos);
        boxplot([tau.All.D1,tau.All.D2]);
        hold on;
        errorbar([1:numel(ROIs)*2], [tau.All.D1_mean, tau.All.D2_mean], [tau.All.D1_se, tau.All.D2_se], '.');
        title('Tau_A for Dimension 1 and Dimension 2 Rank RDM');
        xticklabels([ROIs, ROIs]);
        hold off;
        if svoption
            saveas(fig(i), fullfile(svpath, [svfile, '_D1D2']), 'epsc');
        end
    end
    % Save figs
    if svoption
        savefig(fig, fullfile(svpath, [svfile, surfix, '.fig']));
        close (fig); 
    end
    fprintf('\n Done \n');

    %% Histogram
    if histon && ~maskoption
        for roiIdx=1:numel(ROIs)
            clear gmean ranktauperm
            Tauhist(roiIdx)=figure;
            hist(tau.(ROIs{roiIdx}).tauRaw.permE(:),100);
            hold on;
            scatter(tau.(ROIs{roiIdx}).tauRaw.Ec0, zeros(1,21),'.');
            gmean=mean(mean(tau.(ROIs{roiIdx}).tauRaw.permE'));
            scatter(gmean,0,'x');
            xticks([gmean-0.3, gmean-0.2, gmean-0.1, gmean, gmean+0.1, gmean+0.2, gmean+0.3, gmean+0.4, gmean+0.5]);
            xticklabels({'-0.3', '-0.2', '-0.1', '0.0', '0.1', '0.2', '0.3', '0.4', '0.5'});
            ranktauperm=sort(tau.(ROIs{roiIdx}).tauRaw.permE(:));
            plot([ranktauperm(length(tau.(ROIs{roiIdx}).tauRaw.permE(:))*0.05),ranktauperm(length(tau.(ROIs{roiIdx}).tauRaw.permE(:))*0.95)],[100,100]);
            plot([ranktauperm(length(tau.(ROIs{roiIdx}).tauRaw.permE(:))*0.10),ranktauperm(length(tau.(ROIs{roiIdx}).tauRaw.permE(:))*0.90)],[150,150]);
            plot([ranktauperm(length(tau.(ROIs{roiIdx}).tauRaw.permE(:))*0.01),ranktauperm(length(tau.(ROIs{roiIdx}).tauRaw.permE(:))*0.99)],[50,50]);
            title(ROIs{roiIdx});
            if svoption
                saveas(Tauhist(roiIdx), fullfile(svpath, ['Histogram_Tau_', ROIs{roiIdx}, surfix]), 'epsc'); %surfix
            end
        end
        if svoption
            savefig(Tauhist, fullfile(svpath, ['Histogram_Tau_', surfix, '.fig'])) %surfix
            close(Tauhist);
        end    
    end
end % for model m