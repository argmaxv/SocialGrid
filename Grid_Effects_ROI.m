%% Grid_Effects_ROI
    % visualize grid effects of six-folds symmetry in the brain activity (extracted Beta) from ROIs 

clear; close all; clc;

%% Input 

%%%%%%%%%%%%%%%%%%%%%%%%%%% Update Them %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PhiType={'PhixB', 'PhixD'}; % Phi estimated with cross validation across blocks; Phi estimated with cross validation across days;
typAna0={'All'}; % Analysis including all pairs (both F0-F1 and F0-F2 pairs) 
%typAna0={'All', 'Half', 'Novel', 'F0'};    % Alternatively, one can select other possible
%analyses. See each analysis below.
% All - All inferred vector samples
% Half - Test the grid effects on the bins of 0~pi range only (Control analysis to test whether grid effects were significant only in first 6 bins)
% Novel - Test the grid effects only for the pairs when presented it for the first
% time (therefore it includes lower sample size - some pairs in Day1 scan) (Control test to see whether the grid codes is used for novel inferences)
% F0 - Test the grid effects separately according to its F0 position (Control analysis to test whether the grid effects were selectively employed according to its spatial configurations)
periodrng=[4:8];    % 4~8, range of periodicity
periodicity=6;        % 6, periodicity of interest, six-fold symmetry
svoption=1;           %1, to save the results (in svfolder, below). 0, otherwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[subj, subn] = CallSubj_PS;
PhiROI=ROI.Grid{1};                    % ROI for initial orientation, EC_Grid;
GrdOrientation=PhiROI(1:end-5); % ROI for initial orientation (e.g. 'EC', for save file names)
svfolder=[ProjSet.Phipath, 'GridCategory_', GrdOrientation]; %Where to save the results
zBdatapath=[ProjSet.PhiInfopath, PhiROI, fs, 'zB']; %Where zB (Beta, brain activity) are (generated from EachDec_Beta_Collect.m) 

% Initializing num. of figures.
fig_cat_n=0;
fig_12s_n=0;
fig_12b_n=0;
fig_12s_nv=0;
rec=1; %break time (sec) to make sure figure is drawn
gri=find(ismember(periodrng,periodicity)); %3, position of 6-fold in range of periodicity, 

% Load files generated from text2mat (Information about 'novel pairs' and 'correct decision trials' to make mask for analysis)
load(fullfile([ProjSet.Metapath, PhiROI, fs], 'TrialMask.mat') );
Accuracy=[]; StartingFc=[];
% infomation about trial should be modified in a format of Beta files.
%Beta files includes F1 trial 1 2 ... but also F2 trial 1 2 ... in each
%block. % Therefore, Information about each trial should be repeted once more to cover not only F1 but also F2 time phase
for se=1:info.Nday*info.Nses 
    Accuracy=[Accuracy; repmat(TrMask.accuracynan(TrMask.sess(:,1)==se,:),[2,1])]; % 1 is accurate trial, 0 is inaccurate response.
    StartingFc=[StartingFc; repmat(TrMask.F0(TrMask.sess(:,1)==se,:),[2,1])]; % identity of the F0 of the presenting pair 
end
nvTrial=TrMask.nvTrial; % mask of novel pairs (1 novel, 0 not novel)

%% Main
for Phii=1:numel(PhiType) % Apply differnt Phi (corss validation across blocks and days) to the same brain responses (raw BOLD acitivity, Beta)
    clear AltPeriodicity
	
    % load Theta-Phi computed from GridCategory.m. AltPeriodicity includes
    % the category of bin according the the inferred vector angles, on and
    % off grids (aligned or not to the grid orientation), and periodicities
    % (not only 6-fold but also alternative periodicities).
    load(fullfile('/home/ravic/Pr4.PF/Imaging/Analysis/Analysis Phi/GeneratedModel/EC_Grid/', ['GridCat_', PhiROI, '_', PhiType{Phii}, '.mat']));
    if svoption
        if ~exist(svfolder,'dir')
            mkdir(svfolder);
        end
    end
	if strcmp(PhiType{Phii}, 'PhixB')
        [ROIname, curROI10, indexvoi0]=CallROIzBeta({'EC_xB', 'mPFC_xB','mPFC_GP', 'rTPJ_GP', 'lTPJ_GP'}); % Using Phi xB (across blocks cross validation), test the grid effects of EC, mPFC activity and GP effects encoded in bi-TPJ and mPFC (Growth potential, decision values) 
	elseif strcmp(PhiType{Phii}, 'PhixD')
        [ROIname, curROI10, indexvoi0]=CallROIzBeta({'EC_xD', 'mPFC_xD'}); % Using Phi xD (across days cross validation), test the grid effects of EC, mPFC activity
    end
    clear Onf_com valCtg Crss

    for ri=1:numel(ROIname)
        curROI1=curROI10{ri};
        indexvoi=indexvoi0{ri};
        gridcatfile=dir([zBdatapath, fs, curROI1, '*.mat']);
        load(fullfile(gridcatfile.folder, gridcatfile.name));   % Load the Beta activity extracted from the ROI
        zB=zB.*Accuracy;        % Masking the incorrect trials.
        fprintf('\n ----------- %s ----------- \n', indexvoi);

        for tp=1:numel(typAna0)
        clear typAna
        typAna=typAna0{tp};
        switch typAna              % Define Mask according to the type of analysis
            case 'All'
                clear Themask rep tempmask
                for pr0=1:numel(periodrng)
                    Themask(pr0).onoff=AltPeriodicity(pr0).onoff.*Accuracy;
                end

            case 'Half'
                clear Themask rep tempmask
                for pr0=1:numel(periodrng)
                    tempmask=AltPeriodicity(pr0).onoff.*Accuracy;
                    tempmasknan=nan(size(tempmask));
                    tempmasknan(AltPeriodicity(gri).cat<6)=1; % Category bin 0~5 in six-fold symmetry
                    Themask(pr0).onoff=tempmask.*tempmasknan;
                end

            case 'Novel'
                clear Themask
                for pr0=1:numel(periodrng)
                    clear tempmask
                    tempmask=AltPeriodicity(pr0).onoff.*Accuracy;
                    tempmasknan=nan(size(tempmask));
                    tempmasknan(nvTrial==1)=1;          % Novel trials
                    Themask(pr0).onoff=tempmask.*tempmasknan;
                end
                novelcatmask0=AltPeriodicity(gri).cat.*Accuracy;
                tempmasknan=nan(size(novelcatmask0));
                tempmasknan(nvTrial==1)=1;
                novelcatmask_day1=novelcatmask0.*tempmasknan;

            case'F0'    % Test the grid effects in each of F0 position
                clear Themask tempmask
                for pr0=1:numel(periodrng)
                    Themask(pr0).onoff=AltPeriodicity(pr0).onoff.*Accuracy;
                end
                F0pos=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14]; % Faces shown at the time of F0 (For figure ticklabels)
        end
        
        clear valCtg selDiff
        for sub=1:subn
            if strcmp(typAna,'All')
                for cat=0:2:periodicity*2-1
                        valCtg(cat+1,sub)=nanmean(zB(AltPeriodicity(gri).cat(:,sub)==cat,sub)); %On even
                        valCtg(cat+2,sub)=nanmean(zB(AltPeriodicity(gri).cat(:,sub)==cat+1,sub)); %Off odd
                end
            elseif strcmp(typAna,'Novel')
                for cat=0:2:periodicity*2-1
                        nvalCtg(cat+1,sub)=nanmean(zB(novelcatmask_day1(:,sub)==cat,sub)); %On even
                        nvalCtg(cat+2,sub)=nanmean(zB(novelcatmask_day1(:,sub)==cat+1,sub)); %Off odd
                end                
            else
            end
            
            if strcmp(typAna,'F0')
                clear ontr0 offtr0;
                %pr0=3; %6-fold;
                for F0i=1:numel(F0pos)    
                    ontr0=zB(StartingFc(:,sub)==F0pos(F0i) & Themask(gri).onoff(:,sub)==1,sub); %  On
                    offtr0=zB(StartingFc(:,sub)==F0pos(F0i) & Themask(gri).onoff(:,sub)==0,sub); %  Off
                    F0Pos.(ROIname{ri}).meanOn(sub,F0i)=nanmean(ontr0);                                    %  mean On
                    F0Pos.(ROIname{ri}).meanOff(sub,F0i)=nanmean(offtr0);                                    %	mean Off
                    F0Pos.(ROIname{ri}).selDiff(sub,F0i)=nanmean(ontr0)-nanmean(offtr0);             %  Difference On-Off
                end
            else   
                for pr0=1:numel(periodrng)
                        clear ontr offtr;
                        ontr=zB(Themask(pr0).onoff(:,sub)==1,sub); %On
                        offtr=zB(Themask(pr0).onoff(:,sub)==0,sub); %Off
                        ontr=rmmissing(ontr);
                        offtr=rmmissing(offtr);
                        Onf_com(pr0).(ROIname{ri}).(typAna).meanOff(sub)=nanmean(offtr);
                        Onf_com(pr0).(ROIname{ri}).(typAna).meanOn(sub)=nanmean(ontr);
                        Onf_com(pr0).(ROIname{ri}).(typAna).seOff(sub)=SEM(offtr);
                        Onf_com(pr0).(ROIname{ri}).(typAna).seOn(sub)=SEM(ontr);
                        Onf_com(pr0).(ROIname{ri}).(typAna).selDiff(sub)=nanmean(ontr)-nanmean(offtr);
                end %for pr0=1:numel(periodrng) 
            end  % if periodrng(pr0)==6 && strcmp(typAna,'F0')
        end %for sub
        
        if ~strcmp(typAna,'F0') 
            clear  col
            for pr0=1:numel(periodrng) % 3 indicates Periodicity 6
                clear p stats
                % Pairwise comparisons between alternative periodicity and
                % six-fold ([4 vs. 6] [5 vs. 6] [7 vs. 6]  and [8 vs. 6] )
                Onf_com(pr0).(ROIname{ri}).(typAna).diffMean=nanmean(Onf_com(pr0).(ROIname{ri}).(typAna).selDiff);
                Onf_com(pr0).(ROIname{ri}).(typAna).diffSEM=SEM(Onf_com(pr0).(ROIname{ri}).(typAna).selDiff');
                [~,p,~,stats]=ttest(Onf_com(pr0).(ROIname{ri}).(typAna).selDiff); clear stval; stval=stats.tstat;
                Onf_com(pr0).(ROIname{ri}).(typAna).oneP=p;
                Onf_com(pr0).(ROIname{ri}).(typAna).oneVal=stval;
                for pr1=1:numel(periodrng)
                    clear bin1 bin2 binA binB p stats

                    if pr0==pr1 % exclude comparisons to oneself 
                        Crss.(ROIname{ri}).(typAna).crsP(pr0,pr1)=NaN; 
                        Crss.(ROIname{ri}).(typAna).crsT(pr0,pr1)=NaN;
                    else
                        bin1=Onf_com(pr0).(ROIname{ri}).(typAna).selDiff';
                        bin2=Onf_com(pr1).(ROIname{ri}).(typAna).selDiff';
                        binA=bin1(~isnan(bin1));
                        binB=bin2(~isnan(bin2));                
                        [~,p,~,stats]=ttest(binA,binB);
                        Crss.(ROIname{ri}).(typAna).crsP(pr0,pr1)=p;
                        Crss.(ROIname{ri}).(typAna).crsT(pr0,pr1)=stats.tstat;
                    end %if   
                end %for
                col.means(pr0)=Onf_com(pr0).(ROIname{ri}).(typAna).diffMean;
                col.sems(pr0)=Onf_com(pr0).(ROIname{ri}).(typAna).diffSEM;
                col.pnts(:,pr0)=Onf_com(pr0).(ROIname{ri}).(typAna).selDiff(~isnan(Onf_com(pr0).(ROIname{ri}).(typAna).selDiff))';
                col.prs(:,pr0)=pr0*ones(numel(col.pnts(:,pr0)),1);
            end

          %% On-Off differnt graph
            fig_cat_n=fig_cat_n+1;
            diffBars(fig_cat_n)=figure;
            boxplot(col.pnts,'PlotStyle','compact','OutlierSize',1); %boxplot(col.pnts);
            hold on;
            errorbar([1:1:numel(periodrng)], col.means, col.sems, 'o');
            clear ax
            ax=gca;
            if strcmp(typAna,'All')
                ax.YLim=[-.3, .3];
            else
                ax.YLim=[-1, 1];
            end
            ax.Title.String=[typAna, '   ', indexvoi];
            line(xlim(),[0,0],'LineWidth',1,'Color','k'); 
            pause(rec);
            if svoption
                saveas(diffBars(fig_cat_n), fullfile(svfolder, ['Boxplot_OnOffDiff_',typAna, '_', indexvoi, '_', ROIname{ri}, '.eps']),'epsc');
            end
        end % if ~F0
        
    %% 12 bins
        if strcmp(typAna,'All')
            
            valCtgcol{ri}=valCtg;
            valCtgsel=valCtg;
            fig_12b_n=fig_12b_n+1;
            catbars(fig_12b_n)=figure;
            bar(nanmean(valCtgsel'));
            hold on;
            errorbar([1:12],nanmean(valCtgsel'),SEM(valCtgsel'),'.');
            X_valCtgsel=repmat([1:size(valCtgsel,1)]',[1,size(valCtgsel,2)]);
            scatter(X_valCtgsel(~isnan(valCtgsel)),valCtgsel(~isnan(valCtgsel)),'.');
            clear ax
            ax=gca;
            ax.YLim=[-.3, .3];
            ax.Title.String=[indexvoi];
            pause(rec); 
            fig_12s_n=fig_12s_n+1;
            catboxes(fig_12s_n)=figure;
            boxplot(valCtgsel','PlotStyle','compact','OutlierSize',1); %boxplot(valCtgsel');
            hold on;
            errorbar([1:12],nanmean(valCtgsel'),SEM(valCtgsel'),'o');
            clear ax
            ax=gca;
            ax.YLim=[-.8, .8];
            ax.Title.String=[indexvoi];
            line(xlim(),[0,0],'LineWidth',1,'Color','k');
            pause(rec);
            if svoption
                saveas(catbars(fig_12b_n), fullfile(svfolder, ['Periodicity_Scat_', typAna, '_', indexvoi, '_', ROIname{ri}, '.eps']),'epsc');
                saveas(catboxes(fig_12s_n), fullfile(svfolder, ['Periodicity_Box_' typAna, '_', indexvoi, '_', ROIname{ri}, '.eps']),'epsc');
            end
        end %if strcmp(typAna,'All')
        
        if strcmp(typAna,'Novel')
            nvalCtgcol{ri}=nvalCtg;
            nvalCtgsel=nvalCtg;   
            fig_12s_nv=fig_12s_nv+1;
            catboxes_nv(fig_12s_nv)=figure;
            boxplot(nvalCtgsel','PlotStyle','compact','OutlierSize',1); %boxplot(nvalCtgsel');
            hold on;
            errorbar([1:12],nanmean(nvalCtgsel'),SEM(nvalCtgsel'),'o');
            clear ax
            ax=gca;
            ax.YLim=[-1.5, 1.5];
            ax.Title.String=[indexvoi, ' Novel'];
            line(xlim(),[0,0],'LineWidth',1,'Color','k');
            pause(rec);
            if svoption
                saveas(catbars(fig_12b_n), fullfile(svfolder, ['Periodicity_Scat_', typAna, '_', indexvoi, '_', ROIname{ri}, '.eps']),'epsc');
                saveas(catboxes(fig_12s_n), fullfile(svfolder, ['Periodicity_Box_' typAna, '_', indexvoi, '_', ROIname{ri}, '.eps']),'epsc');
            end
        end
        
        if ~strcmp(typAna,'F0')
            fprintf('%s one sample p-val: %1.4f    T-val: %1.4f   \n' ,typAna, Onf_com(3).(ROIname{ri}).(typAna).oneP, Onf_com(3).(ROIname{ri}).(typAna).oneVal)
        end
	end %for TypAna
    clear typAna_sel
    if ismember('F0', typAna0) 
        typAna_sel=[];
        temp_F0=strfind(typAna0, 'F0');
        for ty0=1:numel(typAna0)
            if isempty(temp_F0{ty0})
                typAna_sel=[typAna_sel, typAna0(ty0)];
            end
        end
    else
        typAna_sel=typAna0;
    end
    
    for t0=1:numel(typAna_sel)        
         fprintf('%s\t Cross p-val:\t %1.4f   %1.4f    %1.4f   %1.4f  \n' ,typAna_sel{t0}, Crss.(ROIname{ri}).(typAna_sel{t0}).crsP(3,1),Crss.(ROIname{ri}).(typAna_sel{t0}).crsP(3,2),Crss.(ROIname{ri}).(typAna_sel{t0}).crsP(3,4),Crss.(ROIname{ri}).(typAna_sel{t0}).crsP(3,5))
         fprintf('%s\t Cross T-val:\t %1.4f   %1.4f    %1.4f   %1.4f  \n' ,typAna_sel{t0}, Crss.(ROIname{ri}).(typAna_sel{t0}).crsT(3,1),Crss.(ROIname{ri}).(typAna_sel{t0}).crsT(3,2),Crss.(ROIname{ri}).(typAna_sel{t0}).crsT(3,4),Crss.(ROIname{ri}).(typAna_sel{t0}).crsT(3,5))
    end
    %% in case 'F0', draw mean on-, off-grid, and difference activity in each F0
    if strcmp(typAna,'F0')
        F0fig(ri)=figure;
        YMean=[nanmean(F0Pos.(ROIname{ri}).meanOn); nanmean(F0Pos.(ROIname{ri}).meanOff); nanmean(F0Pos.(ROIname{ri}).selDiff)]';
        Ysem=[SEM(F0Pos.(ROIname{ri}).meanOn); SEM(F0Pos.(ROIname{ri}).meanOff); SEM(F0Pos.(ROIname{ri}).selDiff)]';
        bar(YMean); hold on;
        x=MultiErrorbar(YMean);
        errorbar(x,YMean, Ysem, 'k', 'linestyle', 'none');
        for f0i=1: numel(F0pos)
            if ~isnan(nanmean(F0Pos.(ROIname{ri}).selDiff(:,f0i)))
                clear p_f0;
                p_f0=rsa.stat.signrank_onesided(F0Pos.(ROIname{ri}).selDiff(:,f0i));
                if p_f0<0.01
                    text(f0i, mean(F0Pos.(ROIname{ri}).selDiff(:, f0i)), '**');
                elseif p_f0<0.05
                    text(f0i, mean(F0Pos.(ROIname{ri}).selDiff(:, f0i)), '*');
                else
                end
            end
        end

        xticklabels(F0pos)
        xlabel('F0');
        ylabel('mean activity');
        title(ROIname{ri});
        pause(rec);
        if svoption
            saveas(F0fig(ri), fullfile(svfolder, ['F0_OnOff_',indexvoi, '_', typAna, '.eps']),'epsc');
        end
    end
end
end
%% save figures in *.fig format and results in *.mat file.
if svoption
    if sum(ismember([{'All'},{'Half'},{'Novel'}], typAna0))
    save( fullfile(svfolder, ['GridCat_Stat_', GrdOrientation, '.mat']), 'Onf_com', 'Crss');
    savefig(diffBars, fullfile(svfolder, ['GridCat_onoff_', GrdOrientation, '.fig']));
    savefig(catbars, fullfile(svfolder, ['GridCat_BarAll_', GrdOrientation, '.fig']));    
    savefig(catboxes, fullfile(svfolder, ['GridCat_BoxAll_', GrdOrientation, '.fig']));  
    savefig(catboxes_nv, fullfile(svfolder, ['GridCat_BoxNv_', GrdOrientation, '.fig']));  
    end
    if sum(ismember({'F0'}, typAna0))
        savefig( F0fig, fullfile(svfolder, ['F0_OnOff_', GrdOrientation, '.fig']) );
    end
    close all;
end

