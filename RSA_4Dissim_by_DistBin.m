% RSA_4 Mean (Dis)similarity by distance bin
% 27 . Aug. 2018
% S.A.Park
% Computing mean representational (dis)similarity across behavioral RDM
% Using normalized [0-1] brain RDM

clear all
close all
clc

%% Setting

[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[subj, nsubj] = CallSubj_PS;
ROIs=ROI.Ana; %ROI list - HCr,HCl, ECr,ECl,M1r, and M1l

%%%%%%%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%%%%% 
models={'Mtv_NMatch14', 'Mtv_AllFaces14'}; % Names of 1st level analysis path (See 'design_name' in Multivariate_Level1.m)
%'Mtv_NMatch14' - 14 faces (except for 1 and 16) while matching the number of presentation of each face
%'Mtv_AllFaces14' - All presentations of 14 faces (except for 1 and 16) including the F0,F1, anf F2 events
rng=2:15; % 2~15 individuals among 16 individuals in social hierarchy
figsave=1; %1 if want to save figures 0 otherwise
normalizeon=1; %Dissimilarity is normalized otherwise 0 for Crossnobis distance (Cross blocks Mahalnobis dist.)
type={'E', 'D1', 'D2'}; %model RDM type
figtype='box'; %'bar';%    %Fig style, Do you want to draw [Boxplot + mean(se)] or [Barplot + scatters] ?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% EucRDM_factorize generates the model RDM (2-D, 2 2-D)
%modelRDMs_cell{1} is E_rdm
%modelRDMs_cell{2} is D1_rdm
%modelRDMs_cell{3} is D2_rdm
[modelRDMs_cell{1}.RDM, modelRDMs_cell{2}.RDM, modelRDMs_cell{3}.RDM]=EucRDM_factorize(rng);

%% Main
for m=1:numel(models)
    fprintf('Drawing %s ...\n', models{m});
    RDMpath=[ProjSet.Respath, models{m}, fs];
    load(fullfile(RDMpath, fname.Darwroi)); % Load RDM normalized from RSA_3TauBar.m
    
    % Set save path
    if figsave
        if strcmp(figtype, 'bar')
            svfile='Bar'; %save file name
        elseif strcmp(figtype, 'box')
            svfile='Box'; %save file name
        end
        svpath=[RDMpath, 'RSA', fs, svfile]; %save path name
        if ~exist(svpath,'dir')
            mkdir(svpath);
        end
    end

    for r=1:numel(ROIs)
        clear subjectRDMs

        if normalizeon
            subjectRDMs=RDMcross.(ROIs{r}).norm;	% Normalized to 0 to 1
            svfix='_Normalized';
        else
            subjectRDMs=RDMcross.(ROIs{r}).mean;   % Crossnobis distance
            svfix='_Mahalanobis';
        end

        for mx=1:numel(modelRDMs_cell)
            clear Temp0 trilMtx Temp x y

    %% Compute mean dissimilarity (y) per bin (x)
            Temp0=modelRDMs_cell{mx}.RDM;
            Temp0=triu(Temp0,1); %select upper triangle matrix
            trilMtx=ones(size(Temp0,1));
            trilMtx=-1*tril(trilMtx,-1);
            Temp=Temp0+trilMtx; %As a reulst, Temp has a behavioral RDM that only include upper triangle mtx and outside of the trinangle is -1  
            xidx=0;
            clear ycnt cntthreshold
            while max(max(Temp))>-1
                xidx=xidx+1; % number of differnt distances in RDM (length of x-axis)
                x(xidx)=max(max(Temp)); %x axis, Distances in Behavioral RDM
                for sidx=1:size(subjectRDMs,3) %Subjects
                    clear BTemp;
                    BTemp0=subjectRDMs(:,:,sidx);
                    BTemp=(triu(BTemp0)+tril(BTemp0)')/2;

                    y(sidx,xidx)=mean(BTemp(Temp==x(xidx) )); %y axis, Dissimilarity when x=x(xidx)
                    ycnt(sidx,xidx)=sum(sum(Temp==x(xidx)));
                end %for subj sidx
                Temp=round(Temp,3);
                curmax=round(x(xidx),3);
                Temp(Temp==curmax)=-1;
            end % while loop

    %% Figure
            cntthreshold=[mean(ycnt)>1];
            cntthreshold(end)=0;
            if mx>1 %Figure window size
                xlength=.06;
                distance2=flip(round(x,0)); % x axis labels 1D
                x_label='1D Distance';
            else
                xlength=.12;
                distance2=flip(round(x.^2,0)); % x axis labels 2D
                x_label='2D Distance^2';
            end
            pos = [((mx-1)*.12) ((r-1)*.15) xlength .15];
            Bars(mx+(r-1)*numel(type))=figure;
            set(gcf, 'Color', ones(1,3), 'Units', 'Normalized', 'Position', pos);
            clear ymean ysem

    %For Crossnobis distance (normalization off)
            clear newY szY newX
            newY=flip(y(:,cntthreshold),2);
            szY=size(newY);
            newX=repmat([1:szY(2)],szY(1),1);

            ymean=mean(newY);
            ysem=SEM(newY);
            if strcmp(figtype, 'box')
                boxplot(flip( y(:,cntthreshold),2 )); hold on;
                errorbar([1:szY(2)], ymean, ysem,'o'); pause(0.5);
            end

            if strcmp(figtype, 'bar')
                bar([1:szY(2)], ymean); hold on;
                errorbar([1:szY(2)], ymean, ysem,'o'); pause(0.5);
                scatter(newX(:), newY(:),'.');
            end

            for d=2:numel(distance2)
                if cntthreshold(d)
                    distance2str{d-1}=num2str(distance2(d));
                end
            end
            xticklabels(distance2str);
            xlabel(x_label);
    
	% Normalization on
            if normalizeon 
                ylabel('Dissimilarity (a.u)'); 
                %ylim([0.8 1.0]);
                if strcmp(figtype, 'box')
                    ylim('auto');
                elseif strcmp(figtype, 'bar')
                    ylim([.8 1]);
                end
            else
                
	%Ylim for Crossnobis distance (Normalization off)
                ylabel('Mahalanobis Dist.');
                switch r
                    case 1
                        %ylim([30.6 31.6]);
                        ylim([30.5 32]);
                    case 2
                        %ylim([32 33]);
                        ylim([31.3 33.7]);
                    case 3
                        %ylim([17 18]);
                        ylim([17 18.5]);
                    case 4
                        %ylim([14.3 15.3]);
                        ylim([14.3 15.7]);
                    case 5
                        %ylim([34.5 35.5]);
                        ylim([34 36.1]);
                    case 6
                        %ylim([34.5 35.5]);
                        ylim([34 36.1]);
                end
            end
            
            clear distance2str
            title([ROIs{r}, '   ', type{mx}]);
            if figsave
                saveas(Bars(mx+(r-1)*numel(type)), fullfile(svpath, [svfile, '_', ROIs{r}, '_MahScat_', type{mx}, svfix]), 'epsc');
            end
        end % for RDM type mx
    end% for ROIs r

    %% Save figures
    if figsave
        if strcat(figtype, 'box')
            savefig(Bars, fullfile(svpath, [svfile, '_MahDis_Box', svfix, '.fig']));
        end
        if strcat(figtype, 'bar')
            savefig(Bars, fullfile(svpath, [svfile, '_MahDis_Scat', svfix, '.fig']));
        end
        close (Bars); 
        fprintf('%s is done!\n', models{m});
    end
end % for models m