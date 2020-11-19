%%Log Data Extraction
%
% by SPark 2017.12.20
%
% Reading *.log and *.txt files from output of Presentation and Creating onset files to set GLM in fMRI analysis
% Please add path 'func' that contains costomized functions before
% execution
%


%% Path, initial variables
clear
close all
clc
[ProjSet, fs, info, ROI, fname]=Call_default_PS;

% % % % % % % % % % % %  ApplyPhi and roitype need to be defined according to which model to test  % % % % % % % % % % % % 

ApplyPhi = 0; % 0, No Phi; 1, Cross blocks Phi acquired in the same day; 2, Cross Days Phi acquired in the differnt day 
roitype = ROI.Grid{1}; %{1}, EC_Grid; {2}, mPFC_Grid; if ApplyPhi==0 then you don't need to define roityp.
svoption=1; % 1 to save mat files, 0 otherwise

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

path=ProjSet.ONSETpath; 
addpath(path);
datafolder=ProjSet.Logpath;
Phipath0=ProjSet.Phipath;   
savepath= ProjSet.Metapath; 
Phipath=ProjSet.PhiInfopath;

if svoption && ~exist(path,'dir') % create directory to save the results
    mkdir(path);
    mkdir(savepath);
    mkdir(Regpath);
    mkdir(Phipath);
end

tunit=info.tunit; % time unit
nblocks=info.Nses; %3 blocks
rng(info.rngset); % random seed for (number matched) downsampling process 
subjlist0 =CallSubj_PS;

%% Find log file
% log files - onset of pulse (to synchronize the fMRI and behavioral data), onset of button response (to take off the motor effects)
% txt files - behavioral data (stimuli, resonses, reaction times etc...)

i=1;
for dayx=1:info.Nday
    if dayx==1
        for sblx=1:numel(subjlist0) 
            subjlist{i}=[info.prefix.day1, subjlist0{sblx}];
            i=i+1;
        end
    else
        for sblx=1:numel(subjlist0)
            subjlist{i}=[info.prefix.day2, subjlist0{sblx}];
            i=i+1;
        end
    end
end

for s=1:numel(subjlist)
    flist(s)=dir( fullfile(ProjSet.Logpath, [subjlist{s}, '*-PS_fMRI_D*.log']) );
end
load(fullfile(path, 'OverlapMat.mat')); % 16x16 matrix which contains the Infomation about which paris has been shown during behavioral training (For the control analysis regarding 'Novel pairs')

%% Get the values from the results of Get-Phi.m

disp(['Phi type is: ' num2str(ApplyPhi)]);
Phi=nan(length(flist), nblocks); % default Phi
switch ApplyPhi
    case 0
        descpt='NoPhi';
        for s=1:length(flist)
            for bl=1:nblocks
                Phi(s,bl)=0; %    Phi=[0,0,0; .... 0,0,0;0,0,0];
            end
        end
        
    case 1
        descpt='PhixB'; %Cross Blocks
        phifile=dir( fullfile([Phipath, roitype, fs], ['*_', roitype, '*.mat']) );
        if numel(phifile)~=1
            error(['There are no or more than 1 Phi data in ', Phipath]);
        end
        load( fullfile([Phipath, roitype, fs], phifile.name) ); % read Phi date generated from Get-Phi function
        fprintf(['\n ', descpt, '\n Read: ', phifile.name, '\n']);
        for s=1:length(flist)   %length(flist): num log files (N days * N Subj)
            clear preSes;
            if strcmp(info.prefix.day1, flist(s).name(1)) %Day 1 scan
                preSes=1;
            elseif strcmp(info.prefix.day2, flist(s).name(1)) % Day 2 scan
                preSes=2;
            else
            error('Check the participants ID.');
            end
            
            % Find correct phi and assign them to Phi(s,bl) where
            % s indicates which data (N days * N Subj)
            % bl indicates which block 
            for bl=1:nblocks
                Phi(s,bl)=Grid6_F01F02_5s(s-(preSes-1)*(length(flist)/2)).vF01F02.ses((preSes-1)*nblocks+bl).Phi_xBCirMean; 
            end %for block
            
            if sum(isnan(Phi(s,:)))
                error(['Beta_Grid6 for ', subjlist{s}, ' is missing']);
            end
        end %for num log files (N days * N Subj)
    
    case 2
        descpt='PhixD'; %Cross Days
        phifile=dir( fullfile([Phipath, roitype, fs], ['*_', roitype, '*.mat']) );
        if numel(phifile)~=1
            error(['There are no or more than 1 Phi data in ', Phipath]);
        end
        load( fullfile([Phipath, roitype, fs], phifile.name) ); % read Phi date generated from Get-Phi function
        fprintf(['\n ', descpt, '\n Read: ', phifile.name, '\n']);
        for s=1:length(flist)   %length(flist): num log files (N days * N Subj)
            clear preSes;
            if strcmp(info.prefix.day1, flist(s).name(1)) %Day 1 scan
                preSes=1;
            elseif strcmp(info.prefix.day2, flist(s).name(1)) % Day 2 scan
                preSes=2;
            else
            error('Check the participants ID.');
            end
            % Find correct phi and assign them to Phi(s,bl) where
            % s indicates which data (N days * N Subj)
            % bl indicates which block 
            Phi(s,1:nblocks)=Grid6_F01F02_5s(s-(preSes-1)*(length(flist)/2)).vF01F02.ses((preSes-1)*nblocks+1).Phi_xDCirMean;
            
            if sum(isnan(Phi(s,:)))
                error(['Beta_Grid6 for ', subjlist{s}, ' is missing']);
            end
        end %for num log files (N days * N Subj)
        
end %switch

if sum(~isnan(Phi(:)))~=length(flist)*nblocks
    msg = 'Error: %% Num of Phi doent match with Num of subjects %%';
    error(msg)
end

%% log and txt files --> mat data 
nvCheckMat=[];
for s=1:length(flist) % number of log files (N days x N subjects)

%% Reading .log files and get information about the onset times of pulses, block, and button presses   
% Converting log files to mat    
    subj=flist(s).name(1:4); %participant ID 
    sublist{s,1}=subj;
    
    disp(['Process the log file: ' flist(s).name]);
    logstr=importdata( fullfile(datafolder, flist(s).name) );
    logcell=logstr.textdata(4:end,1:5);
    log(:,1)=str2double(logcell(:,2));% Event Index
    log(:,2)=str2double(logcell(:,4));% Event Type Code (e.g. pulse)
    log(:,3)=str2double(logcell(:,5));% Event Onset (time)
    
    % Initial variables for log file reading
    brk=0;
    nrsp=0;
    flag=1;
    starttime=0;
    
    for t=1: length(log)
        if log(t,2)==111 && flag==1% '111' is the specific code to inform the starting time of each block.
            brk=brk+1;
            starttime=log(t,3);
            flag=0;
            
        elseif starttime~=0 && strcmp(logcell(t,3),'Response') && log(t,2)~=3 %when participants press button after scanning but not 3 (3 is pressed by the experimentor)
            nrsp=nrsp+1; %num response
            se(brk).Btn.ons(nrsp,1)=(log(t,3)-starttime)/(10*tunit); % response onset
            se(brk).Btn.typ(nrsp,1)=(log(t,2)-1); % which button 0:left; 1:right
            
        elseif log(t,2)==3 %When the experimenter launched the new block after short break
            flag=1;
            nrsp=0;
            
        else
        end
    end
    clear log;
    clear logcell;
    
%% Initiating variables to save the vector distances, vector angles, and decision values accoring to in which space (coordinate) they are computed.     
    % Initiating variables when reading a new log among [N Days x N Subjects] files 

    tye={'Euc', 'Plc', 'Sig', 'Mds'};  %Euclidean; Placement data; Sigmoid; MDS
    T=initRegs(tye); %T.Euc.GP1=[]; % GP_F01 ...

%% Checking novel pairs: To see whether the pair was presented before the current trial or not

    nvCnt.(flist(s).name(1:4))=0;
    nvTrial.(subj).xse=[];
    clear novelpairs
    if ~isfield(nvCheckMat,flist(s).name(1:4)) % Matrix for checking novel pairs: To see whether the pair was presented before the current trial or not
        novelpairs=zeros(info.Nface,info.Nface);
        countpairs=zeros(info.Nface,info.Nface);
        nvCheckMat.(flist(s).name(1:4))=zeros(info.Nface,info.Nface);
    else
        novelpairs=nvCheckMat.(flist(s).name(1:4));
        countpairs=nvCountMat.(flist(s).name(1:4));
    end
    
%% Converting response data (txt files) to mat
    % Extracting infomation from txt files to set up the GLM
    
    flist1=dir( fullfile(datafolder, ['f.', subj '*.txt']) ); % Behavioral data txt files
    for ses=1:length(flist1)


        fnametxt=flist1(ses).name;
        disp(['Processing txt file:' fnametxt]);
        data=importdata( fullfile(datafolder,fnametxt) ); 
        clear curses; %Data of the current session
        all_curses=data.data;   %Focusing on the current session. Data of each session was recorded in differnt txt files.
        correct_curses=data.data(data.data(:,13)==1,:); %Trials in which participants made correct decisions
        clear curses;
        curses=all_curses;
        condition={'F0a',	'F1',	'F0b',	'F2',   'Fons', 	'F0aOns',   'F1place',  'F2place',  'F1sig' ,'F2sig',   'F1mds',    'F2mds'}; typ=1;
%% Onset time
        %  txt files order
        % 1. Stime	2. Ses	3. Trial	4. DcOnset	5. RT	6. Face0	7. Face1	8. Face2	9. F01OS	10. F1OS	11. F02OS	12. Dec	13. Acc	14. Difficulty	15. GridType	16. HeurisDif	17. Expe
        F1onscnt=10;
        F2onscnt=4; 
        Accurate=13;
        se(ses).(condition{typ,1}).ons=curses(:,9)/tunit; % F0 1st onset 
        se(ses).(condition{typ,3}).ons=curses(:,11)/tunit; %  F0 2nd onset
        se(ses).(condition{typ,2}).ons=curses(:,F1onscnt)/tunit; %  F1 onset
        se(ses).(condition{typ,4}).ons=curses(:,F2onscnt)/tunit; % F2 onset
        se(ses).(condition{typ,4}).dur=curses(:,5)/tunit; % F2 duration (RT)
        %correct trials
        se(ses).(condition{typ,2}).correct=curses(:,Accurate);
        se(ses).(condition{typ,4}).correct=curses(:,Accurate);

        clear day reg_now reg0_now
        if strcmp(info.prefix.day1, fnametxt(3))
            day=0;
        elseif strcmp(info.prefix.day2, fnametxt(3)) 
            day=1;
        else
            error('Check txt file name.');
        end

%% Rank, Grid, and Probability
	% position based on Euclidean space
        F0rank=faceinfo(curses(:,6));
        F0C=F0rank(:,1);    F0C_Euc=F0rank(:,1);
        F0P=F0rank(:,2);    F0P_Euc=F0rank(:,2);
        F1rank=faceinfo(curses(:,7));
        F1C=F1rank(:,1);    F1C_Euc=F1rank(:,1);
        F1P=F1rank(:,2);    F1P_Euc=F1rank(:,2);
        F2rank=faceinfo(curses(:,8));
        F2C=F2rank(:,1);    F2C_Euc=F2rank(:,1);
        F2P=F2rank(:,2);    F2P_Euc=F2rank(:,2);

    %position based on Placement
        F0rank_Plc=faceplace(curses(:,6), 'p', str2double(fnametxt(5:6)));
        F0P_Plc=F0rank_Plc(:,1);F0C_Plc=F0rank_Plc(:,2);                
        F1rank_Plc=faceplace(curses(:,7), 'p', str2double(fnametxt(5:6)));                
        F1P_Plc=F1rank_Plc(:,1);F1C_Plc=F1rank_Plc(:,2);
        F2rank_Plc=faceplace(curses(:,8), 'p', str2double(fnametxt(5:6)));
        F2P_Plc=F2rank_Plc(:,1);F2C_Plc=F2rank_Plc(:,2);

    %position based on a Sigmoid curve
        F0rank_Sig=faceplace(curses(:,6), 's');
        F0P_Sig=F0rank_Sig(:,1);F0C_Sig=F0rank_Sig(:,2);                
        F1rank_Sig=faceplace(curses(:,7), 's');                
        F1P_Sig=F1rank_Sig(:,1);F1C_Sig=F1rank_Sig(:,2);
        F2rank_Sig=faceplace(curses(:,8), 's');
        F2P_Sig=F2rank_Sig(:,1);F2C_Sig=F2rank_Sig(:,2);

    %position based on MDS
        F0rank_Mds=faceplace(curses(:,6), 'm', str2double(fnametxt(5:6)));
        F0P_Mds=F0rank_Mds(:,1);F0C_Mds=F0rank_Mds(:,2);                
        F1rank_Mds=faceplace(curses(:,7), 'm', str2double(fnametxt(5:6)));                
        F1P_Mds=F1rank_Mds(:,1);F1C_Mds=F1rank_Mds(:,2);
        F2rank_Mds=faceplace(curses(:,8), 'm', str2double(fnametxt(5:6)));
        F2P_Mds=F2rank_Mds(:,1);F2C_Mds=F2rank_Mds(:,2);

	% Loading the pairs presented for the behavioral training of each participant to check which pairs were novel.
        clear  overlap_pairs
        overlap_pairs=overlappairs.(subj(2:end));

        for nv=1:size(curses,1)
            % novelpairs has boolean value: 1, Novel; 0, Non-novel
            if novelpairs(curses(nv,6),curses(nv,7))==0 && overlap_pairs(curses(nv,6),curses(nv,7))==0
                novelpairs(curses(nv,6),curses(nv,7))=1;
                nvTrial.(subj).se(ses).F1(nv)=1;
            else
                nvTrial.(subj).se(ses).F1(nv)=0;
            end

            if novelpairs(curses(nv,6),curses(nv,8))==0 && overlap_pairs(curses(nv,6),curses(nv,8))==0
                novelpairs(curses(nv,6),curses(nv,8))=1;
                nvTrial.(subj).se(ses).F2(nv)=1;
            else
                nvTrial.(subj).se(ses).F2(nv)=0;
            end           
            % countpairs contains the numbers of presentation of each pairs
            countpairs(curses(nv,6),curses(nv,7))=countpairs(curses(nv,6),curses(nv,7))+1;
            countpairs(curses(nv,6),curses(nv,8))=countpairs(curses(nv,6),curses(nv,8))+1;
        end

    % myGrid returns radian and Euc. dist
    % Therefore, myGrid(c-a,d-b): the vector from F(a,b) to F(c,d)
        [GrF01, EcF01]=myGrid(F1C-F0C, F1P-F0P);  % F1 -> F0 
        [GrF02, EcF02]=myGrid(F2C-F0C, F2P-F0P);  % F2 -> F0        
        [GrF12, EcF12]=myGrid(F2C-F1C, F2P-F1P);  % F2 -> F1
        [GrD, EcD]=myGrid(F2C-F0C+F1C-F0C, F2P-F0P+F1P-F0P); % Decision vec. F2->F0 - F1->F0  

        for t=1:numel(tye) 
            eval(['[GrF01_', tye{t}, ', EcF01_', tye{t}, ']=myGrid(F1C_', tye{t}, '-F0C_', tye{t}, ', F1P_', tye{t}, '-F0P_', tye{t}, ');'])  % F1 -> F0
            eval(['[GrF02_', tye{t}, ', EcF02_', tye{t}, ']=myGrid(F2C_', tye{t}, '-F0C_', tye{t}, ', F2P_', tye{t}, '-F0P_', tye{t}, ');'])  % F2 -> F0        
            eval(['[GrF12_', tye{t}, ', EcF12_', tye{t}, ']=myGrid(F2C_', tye{t}, '-F1C_', tye{t}, ', F2P_', tye{t}, '-F1P_', tye{t}, ');'])  % F2 -> F1
        end

    % Initiating the variables
        clear GP1 GP2 GrF1GP GrF2GP EcF1GP EcF2GP;
        clear GP1_Euc GP2_Euc GrF1GP_Euc GrF2GP_Euc EcF1GP_Euc EcF2GP_Euc;
        clear GP1_Plc GP2_Plc GrF1GP_Plc GrF2GP_Plc EcF1GP_Plc EcF2GP_Plc;
        clear GP1_Sig GP2_Sig GrF1GP_Sig GrF2GP_Sig EcF1GP_Sig EcF2GP_Sig;
        clear GP1_Mds GP2_Mds GrF1GP_Mds GrF2GP_Mds EcF1GP_Mds EcF2GP_Mds;

    % Computing the decision variable (growth potential, GP) and the Euclidean
    % distance of GP vectors
        for idx=1: length(curses)
            GP1(idx)=max(F0C(idx), F1C(idx))*max(F0P(idx), F1P(idx));
            GP2(idx)=max(F0C(idx), F2C(idx))*max(F0P(idx), F2P(idx));

            for t=1:numel(tye)
                eval(['GP1_', tye{t},'(idx)=max(F0C_', tye{t},'(idx), F1C_', tye{t},'(idx))*max(F0P_', tye{t},'(idx), F1P_', tye{t},'(idx));']);
                eval(['GP2_', tye{t},'(idx)=max(F0C_', tye{t},'(idx), F2C_', tye{t},'(idx))*max(F0P_', tye{t},'(idx), F2P_', tye{t},'(idx));']);
            end
        end

        for t=1:numel(tye)
            eval([' T.(tye{t}).GP1=[T.(tye{t}).GP1; transpose(GP1_', tye{t},')]; ']); % GP_F01
            eval([' T.(tye{t}).GP2=[T.(tye{t}).GP2; transpose(GP2_', tye{t},')]; ']); % GP_F02
            eval([' T.(tye{t}).GrF01=[T.(tye{t}).GrF01; transpose(GrF01_', tye{t},')]; ']); % Gr_F01
            eval([' T.(tye{t}).GrF02=[T.(tye{t}).GrF02; transpose(GrF02_', tye{t},')]; ']); % Gr_F02
            eval([' T.(tye{t}).EcF01=[T.(tye{t}).EcF01; transpose(EcF01_', tye{t},')]; ']); % Ec_F01
            eval([' T.(tye{t}).EcF02=[T.(tye{t}).EcF02; transpose(EcF02_', tye{t},')]; ']); % Ec_F02
        end

    % Compute the probability of winning of F0F1 or F0F2 based on GP_F0F1 and F0 position
        % how much the decision is predictable before presenting F2.
        % The number of faces of F2 that makes F01 wins or F02 wins. 
        % normalized from 0-1.
        
        clear pF1win pF2win uF1win;
        for idx=1: length(curses) % for 1
            pF1win(idx)=0; pF2win(idx)=0; %probability of winning
            for ic=1:4 % for 2
                for ip=1:4 % for 3
                    if (F0C(idx)==ic) && (F0P(idx)== ip) % if 1
                    else
                        if GP1(idx)<max(F0C(idx), ic)*max(F0P(idx), ip) %if 2
                            pF2win(idx)=pF2win(idx)+1;
                        elseif GP1(idx)>max(F0C(idx), ic)*max(F0P(idx), ip)
                            pF1win(idx)=pF1win(idx)+1;
                        else
                        end %if 2
                    end% if 1
                end % for 3
            end % for 2
        end % for 1
            
        % Uncertainty : 1 indicates highly uncertain: p(F1 win)=p(F2 win), 0 indicates it is certain - 1/exp(1) (~0.367) indicates no need to see F2
        % Because that p(GPF0F1)=p(GPF0F2), 1 = p(GPF0F1>GPF0F2) + p(GPF0F1<GPF0F2) + p(GPF0F1=GPF0F2);
        % Uncertainty = 1/ (p(GPF0F1>GPF0F2)-p(GPF0F1<GPF0F2));
        pF1win=pF1win/(16-2); %14 faces except for f0 anf f1
        pF2win=pF2win/(16-2); %14 faces except for f0 anf f1
        uF1win=1/(1-(1/exp(1))).*(1./exp(abs(pF1win-pF2win)) - 1/exp(1));
            
%% Getting onsets of each of all faces presentation for multivariate analysis
%% 1 Onset of every face presentation (Category: Face1, Face2, ... Face16)       
        clear EachFace Img F0a Faces F0aOns F0Ons FaceNMatch lots F0
        Img=[];
        Img(:,2)=[curses(:,9)/tunit; curses(:,10)/tunit;curses(:,11)/tunit;curses(:,4)/tunit]; % onset times
        Img(:,1)=[curses(:,6); curses(:,7);curses(:,6);curses(:,8)]; % what faces   
        Img=sortrows(Img,1);
        F0(:,1)=[curses(:,6);curses(:,6)];
        F0(:,2)=[curses(:,9);curses(:,11)]/tunit;
        F0=sortrows(F0,1);

        EachFace=[];
        Faces={};
        preCount=0;
        for num=1:info.Nface
            count=0;
            for len=1:length(Img(:,1))
               if Img(len,1)==num
                   count=count+1;
                   Faces.(['Face', num2str(num)]).ons(count,1)=Img(len,2); % onset of every face presentation (each face category)
               end
            end
            EachFace.count(num)=count;
            EachFace.start(num)=preCount;
            preCount=preCount+count;
        end
        clear count; clear num; clear len

%% 1-2 Onset of every face presentation of each trial            
        AllFace=sortrows(Img,1);
        for lst=1:length(Img)
            EachFace.name{lst,1}=['F',sprintf('%02d',Img(lst,1)), '_', sprintf('%03d',lst)]; % onset of every face presentation of each trial
            EachFace.Ons(lst,1)=Img(lst,2);
        end

%% 2 Onset of every F0 presentation (Category: Face2, Face3 ... )
        F0aOns={};
        F0Ons={};
        for num=1:info.Nface
            counta=0;
            count0=0;
            for len=1:length(F0(:,1))
                if F0(len,1)==num
                    count0=count0+1;
                    F0Ons.(['Face', num2str(num)]).ons(count0,1)=F0(len,2);  % Onset of F0 presentation 
                end
            end
        end

%% 3. Number of sample size matching (Except for Face 1 and Face 16)
    % If a face is presented more frequent than others, downsample them by
    % selecting random trials of them to match the sample sizes to be the same
    % for all 14 faces
        FaceNMatch=[];
        for fn=2:15 %Except for Face 1 and Face 16
            clear lots
            lots=[randperm(EachFace.count(fn),min(EachFace.count(2:15)))]+EachFace.start(fn); 
            Fmatch.(['Face', num2str(fn)])=EachFace.Ons(lots);
            FaceNMatch=[FaceNMatch;lots'];
        end

        FaceNMatch=sortrows(FaceNMatch);
        EachFace.NmatchOns=EachFace.Ons(FaceNMatch);
        EachFace.NmatchName=EachFace.name(FaceNMatch);
        EachFace.Fmatch=Fmatch;
        clear count; clear num; clear len        

    %% Assign them to parametric modulators (PM) in each phase

%% 1. At the time of F0a (the first F0) presentation
            n=1; se(ses).(condition{typ,1}).pm(:,n)=F0C;                    % 1 competence
            n=n+1; se(ses).(condition{typ,1}).pm(:,n)=F0P;                % 2 popularity

            [se(ses).(condition{typ,1}).CoR, se(ses).F0a.CoP]=corr(se(ses).F0a.pm);
            se(ses).(condition{typ,1}).pmidx={'1 CF0'; '2 PF0';};

%% 2. At the time of F1 presentation (Euclidean space)
            n=1; se(ses).(condition{typ,2}).pm(:,n)=F1C;                           % 1 competence rank
            n=n+1; se(ses).(condition{typ,2}).pm(:,n)=F1P;                       % 2 popularity rank
            n=n+1; se(ses).(condition{typ,2}).pm(:,n)=GP1;                      % 3 decision value, GP01
            n=n+1; se(ses).(condition{typ,2}).pm(:,n)=pF1win;                  % 4 probability F1 win
            n=n+1; se(ses).(condition{typ,2}).pm(:,n)=pF2win;                  % 5 probability F2 win
            n=n+1; se(ses).(condition{typ,2}).pm(:,n)=uF1win;                  % 6 uncertainty given F0 and F1
            n=n+1; se(ses).(condition{typ,2}).pm(:,n)=GrF01;                   % 7 Rad F01 vec. (Theta)
            n=n+1; se(ses).(condition{typ,2}).pm(:,n)=GrF01-Phi(s,ses);   % 8 Rad-Phi F01 vec (Theta - Phi)
            n=n+1; se(ses).(condition{typ,2}).pm(:,n)=EcF01;                    % 9 Euc. dist. F01
            
            [se(ses).(condition{typ,2}).CoR, se(ses).(condition{typ,2}).CoP]=corr(se(ses).(condition{typ,2}).pm); %Correation between potential parametric modulators
            se(ses).(condition{typ,2}).pmidx={'1 CF1'; '2 PF1'; '3 GPF01'; '4 pF1win'; '5 pF2win'; '6 UF1'; '7 R0F01'; '8 RPhiF01'; '9 EucF01'}; %; '10 RiPvF01'; '11 RF1GP'; '12 RPF1GP'; '13 EF1GP';'14 Novel1';};

%% 3. At the time of F1 presentation (Alternative spaces)
            for t=2:numel(tye) %Alternative space t strats from 2 because 1 is Euclidean space.
                cni=[7,9,11]; % according to how the condition was defiend         condition={'F0a',	'F1',	'F0b',	'F2',   'Fons',	'F0aOns',   'F1place',  'F2place',  'F1sig' ,'F2sig',   'F1mds',    'F2mds'};
                eval(['n=1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=F1C_',tye{t}, ';']);                          % 1 competence rank
                eval(['n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=F1P_',tye{t}, ';']);                      % 2 popularity rank
                eval(['n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=GP1_',tye{t}, ';']);                     % 3 decision value GP01
                eval(['n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=GrF01_',tye{t}, '-Phi(s,ses);']); % 4 Rad-Phi F01 vec (Theta - Phi)
                eval(['n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=EcF01_',tye{t}, ';']);                  % 5 Euc. dist. F01 in subjective space
                n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=GP1;                                                   % 6 EucSpace_GP01
                n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=GrF01-Phi(s,ses);                                % 7 EucSpace_Rad-Phi F01 vec (Theta - Phi)
                n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=EcF01;                                                 % 8 EucSpace_Euc dist. F01
                
                [se(ses).(condition{typ,cni(t-1)}).CoR, se(ses).(condition{typ,cni(t-1)}).CoP]=corr(se(ses).(condition{typ,cni(t-1)}).pm); %Correation between potential parametric modulators
                se(ses).(condition{typ,cni(t-1)}).pmidx={'1 CF1','2 PF1','3 GPF01_sub','4 RPhiF01_sub','5 EvF01_sub','6 GPF01_Euc','7 RPhiF01_Euc','8 EvF01_Euc'};
            end

%% 4. At the time of F0b (the second F0) presentation
             n=1; se(ses).(condition{typ,3}).pm(:,n)=F0C;                            % 1 competence rank
             n=n+1; se(ses).(condition{typ,3}).pm(:,n)=F0P;                       % 2 popularity rank
             n=n+1; se(ses).(condition{typ,3}).pm(:,n)=GP1;                      % 3 decision value GP01
             n=n+1; se(ses).(condition{typ,3}).pm(:,n)=pF1win;                  % 4 probability F1 win
             n=n+1; se(ses).(condition{typ,3}).pm(:,n)=pF2win;                  % 5 probability F2 win
             n=n+1; se(ses).(condition{typ,3}).pm(:,n)=uF1win;                  % 6 uncertainty given F0 and F1
             n=n+1; se(ses).(condition{typ,3}).pm(:,n)=GrF01;                   % 7 Rad F01 vec. (Theta)
             n=n+1; se(ses).(condition{typ,3}).pm(:,n)=GrF01-Phi(s,ses);  % 8 Rad-Phi F01 vec (Theta - Phi)
             n=n+1; se(ses).(condition{typ,3}).pm(:,n)=EcF01;                   % 9 Euc. dist. F01

             [se(ses).(condition{typ,3}).CoR, se(ses).(condition{typ,3}).CoP]=corr(se(ses).(condition{typ,3}).pm); %Correation between potential parametric modulators
             se(ses).(condition{typ,3}).pmidx={'1 CF0'; '2 PF0'; '3 GPF01'; '4 pF1win'; '5 pF2win'; '6 UF1'; '7 R0F01'; '8 RPhiF01'; '9 EvF01'}; 

%% 5. At the time of F2 presentation (Euclidean spaces)
            n=1; se(ses).(condition{typ,4}).pm(:,n)=F2C;                            % 1 competence rank
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=F2P;                        % 2 popularity rank
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=GP2;                       % 3 decision value GP02
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=curses(:,14);           % 4 Easiness, abs(GP1-GP2)
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=uF1win;                   % 5 uncertainty by F1
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=GrF02;                    % 6 Rad F02 vec. (Theta)
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=GrF02-Phi(s,ses);    % 7 Rad-Phi F02 vec (Theta - Phi)
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=EcF02;                    % 8 Euc F0F2
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=GrF12;                    % 9 Rad F1->F2
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=GrF12-Phi(s,ses);   % 10 Rad-Phi F12 vec
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=EcF12;                    % 11 Euc F1->F2
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=GrD;                       % 12 Rad Dec vec (vecF0F1 - vecF0F2)
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=GrD-Phi(s,ses);      % 13 Rad-Phi Dec vec
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=EcD;                       % 14 Euc F0F1 - F0F2
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=curses(:,17)-1;       %15  Expected decision
            n=n+1; se(ses).(condition{typ,4}).pm(:,n)=EcF01;                    %16 Euc F0F1

            [se(ses).(condition{typ,4}).CoR, se(ses).(condition{typ,4}).CoP]=corr(se(ses).(condition{typ,4}).pm);%Correation between potential parametric modulators
            se(ses).(condition{typ,4}).pmidx={  '1 CF2'; '2 PF2'; '3 GPF02'; '4 Easiness'; '5 UF1'; '6 R0F02'; '7 RPhiF02'; '8 EucF02';...
                                                                        '9 R0F12'; '10 RPhiF12'; '11 EucF12'; '12 R0Dec'; '13 RPhiDec'; '14 EucDec'; '15 ExpDec'; '16 Euc F0F1'};

%% 6. At the time of F2 presentation (Alternative spaces)
        for t=2:numel(tye)
            cni=[8,10,12]; % according to the definition of conditoin         condition={'F0a',	'F1',	'F0b',	'F2',   'Fons',	'F0aOns',   'F1place',  'F2place',  'F1sig' ,'F2sig',   'F1mds',    'F2mds'};
            eval(['n=1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=F2C_',tye{t}, ';']);                             % 1 competence
            eval(['n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=F2P_',tye{t}, ';']);                         % 2 popularity
            eval(['n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=GP2_',tye{t}, ';']);                        % 3 GP02
            eval(['n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=GrF02_',tye{t}, '-Phi(s,ses);']);    % 4 Rad-Phi F02 vec (Theta - Phi)
            eval(['n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=EcF02_',tye{t}, ';']);                     % 5 Euc dist. F02 in subjective space
            eval(['n=n+1; se(ses).(condition{typ,8}).pm(:,n)=abs(GP2_',tye{t}, '-GP1_',tye{t},');']);  % 6 Easiness, abs(GP1-GP2)
            n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=GP2;                                                       % 7 GP02 _ Euc. space
            n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=GrF02-Phi(s,ses);                                   % 8 Rad-Phi F02 vec (Theta - Phi)_ Euc. space
            n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=EcF02;                                                    % 9 Euc F02_ Euc. space
            n=n+1; se(ses).(condition{typ,cni(t-1)}).pm(:,n)=curses(:,14);                                          %10 Easiness, abs(GP1-GP2)_ Euc. space
            
            [se(ses).(condition{typ,cni(t-1)}).CoR, se(ses).(condition{typ,cni(t-1)}).CoP]=corr(se(ses).(condition{typ,cni(t-1)}).pm);%Correation between potential parametric modulators
            se(ses).(condition{typ,cni(t-1)}).pmidx={'1 CF1','2 PF1','3 GPF02_sub','4 RPhiF01_sub','5 EvF01_sub','6 Easiness_sub','7 GPF01_Euc','8 RPhiF01_Euc','9 EvF01_Euc','10 Easiness_Euc'};
        end
        
%% 7. Onsets of multivariate analyses
        se(ses).([condition{typ,5}])=Faces;
        %se(ses).([condition{typ,6}])=F0aOns;
        se(ses).F0ab=F0Ons; %combining F0a and F0b
        se(ses).EachFace=EachFace; % each trial is modeled separately

%% 8. Meta infomation
        trialInfo.sess(ses).sub(s).accuracy=curses(:,13);
        trialInfo.sess(ses).sub(s).Easiness=curses(:,14);
        trialInfo.sess(ses).sub(s).RT=curses(:,5)/tunit;
        trialInfo.sess(ses).sub(s).F0=curses(:,6);
        trialInfo.sess(ses).sub(s).F1=curses(:,7);
        trialInfo.sess(ses).sub(s).F2=curses(:,8);

        trialInfo.sess(ses).sub(s).F01GP=GP1;
        trialInfo.sess(ses).sub(s).F01Euc=EcF01;
        trialInfo.sess(ses).sub(s).F02GP=GP2;
        trialInfo.sess(ses).sub(s).F02Euc=EcF02;

        trialInfo.sess(ses).sub(s).F01Rad=GrF01'; 
        trialInfo.sess(ses).sub(s).F02Rad=GrF02';
        trialInfo.sess(ses).sub(s).F12Rad=GrF12';
        sub(s).se(ses)=se(ses); % assign to each participant
	end % for session (blocks)
    
%% Creating a mask of novel trials.
     % Fit to the format of the extracted Beta files
     % (48 trials x 2 Phases (F1 and F2) + 7 (gab between sessions)) x
     % nBlock = 309
    for ses=1:length(flist1)
        nvTrial.(subj).xse=[nvTrial.(subj).xse; [nvTrial.(subj).se(ses).F1, nvTrial.(subj).se(ses).F2, nan(1,7)]']; %nan(1,7) is for the gap between sessions (6 motion regressors and button press)
    end
    
%% To what extent GP,cos Angles, and Euclidean distances in subjective spaces correlate with those of the 4x4 Euclidean space.
    for t=2:numel(tye) 
        [placeR.(tye{t}).(flist(s).name(1:4)).r_GP, ~]=corr([T.Euc.GP1;T.Euc.GP2], [T.(tye{t}).GP1;T.(tye{t}).GP2]);
        [placeR.(tye{t}).(flist(s).name(1:4)).r_cosA, ~]=corr(cos([T.Euc.GrF01;T.Euc.GrF02]), cos([T.(tye{t}).GrF01;T.(tye{t}).GrF02]));
        [placeR.(tye{t}).(flist(s).name(1:4)).r_Ec, ~]=corr([T.Euc.EcF01;T.Euc.EcF02], [T.(tye{t}).EcF01;T.(tye{t}).EcF02]);
    end
    
%% Collecting meta infomations (regressors, novel pair mask, trial info)
    subfiledType=fieldnames(T.Euc);     % 'GP1', 'GP2', 'GrF01', 'GrF02', 'EcF01', 'EcF02'
    for t=1:numel(tye)
        for st=1:numel(subfiledType)
                Matric_collect.(tye{t}).(subfiledType{st})(:,s)=T.(tye{t}).(subfiledType{st});
        end
    end    
    nvCheckMat.(flist(s).name(1:4))=novelpairs;                    % novelpairs - Boolean 1, Novel; 0:Non-novel %(flist(s).name(1:4)) is the subject ID
    nvCountMat.(flist(s).name(1:4))=countpairs;                    % countpairs - The number of presentation of each pair

%% Save results    
    if svoption && ApplyPhi==0 
        svfolder=[savepath, fs, 'NoPhi'];
    elseif svoption && (ApplyPhi==1 || ApplyPhi==2)
        svfolder=[savepath, fs, roitype, fs, descpt];
    else %svoption==0
    end

    clear se;
    savename=['Onset_phi_', subj '.mat'];
    Session=sub(s);
    if svoption
        if ~exist(svfolder, 'dir')  
            mkdir(svfolder);
        end
        save([svfolder, fs, savename], 'Session');
    end
end% subject
    
%% Creating a mask of novel trials (combining two days scans)
     % Fit to the format of the extracted Beta files
     % 618 x nSubjects matrix
     % ((48 trials x 2 Phases (F1 and F2) + 7 (gab between sessions)) x
     % nBlock) x 2Days = 618
nvTrial_Collection=[];
for s=1:length(flist)/2
    subj=flist(s).name(1:4);
    nvTrial_Collection=[nvTrial_Collection, [ nvTrial.([info.prefix.day1, flist(s).name(2:4)]).xse; nvTrial.([info.prefix.day2, flist(s).name(2:4)]).xse ]];
end

%% Save metafile. Metafile combines all basic behavioral data information across all blocks, days, and subjects into a matrix
if svoption
    %cd (svfolder);
    savename=fullfile(svfolder, ['ProcessHistory_' descpt, '.mat']);
    save(savename, 'sublist', 'sub', 'ApplyPhi','placeR', 'nvCheckMat', 'nvCnt','nvTrial_Collection', 'trialInfo','Phi','Matric_collect');
    clear all
end


