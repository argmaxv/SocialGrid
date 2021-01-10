function  [BetaInfoPM, BetaPMNum, BetaFilepath]=Get_Beta_Labels(model_name, PM_interest)
% Read SPM.mat from the GLM 1st Lv analysis [model_name]
% and get the index of beta files to extract. To get the beta infomation
% modulated by specific parametric modulator then input [PM_interest] as
% well e.g. Get_BetaInfo_osf('EachDec') or
% Get_BetaInfo_osf('EachDec_GP', 'GP')
% return BetaFilepath - the file path, and BetaInfoPM - Beta file information and
% BetaPMNum - the index of beta file to extract.

    if nargin<2 || isempty(PM_interest)
        PM_interest='bf'; % If there is no specific parametric modulator (PM) then get the information of F1 and F2
    end
    [ProjSet, fs, ~, ~, ~]=Call_default_PS;
    [subj, ~] = CallSubj_PS;
    betapath=[ProjSet.ANA1path, model_name]; % Where Beta File and SPM.mat are located
    
        clear SPM 
        %clear B_Raw B_Sv B_Cv B_GP B_Euc BetaInfo
        sub=4; % Get Beta info only from one subject considering all GLM specifications are the same across participants
        load([betapath, fs, subj{sub}, fs, 'SPM.mat']);
        nBeta=size(SPM.Vbeta, 2); % Number of Beta files

        Block=1;
        for bi=1:nBeta
            clear Bdescrip
            % From Bdescrip in SPM.mat get the information about the Beta file, e.g.  Bdescrip='spm_spm:beta (0001) - Sn(1) F101*bf(1)' 
            Bdescrip=SPM.Vbeta(bi).descrip;

            if bi==1 %After get the infomation from the first Beta file description in Bdescrip, use the infomation
                Bi_pos=strfind(Bdescrip, '0001');   % position of Beta file index
                Bidx=Bdescrip(Bi_pos:Bi_pos+3);  % 000X

                Ev_pos=strfind(Bdescrip, 'F1');       % Position of event
                EventType=Bdescrip(Ev_pos:Ev_pos+1); % F1 or F2
                EventName=Bdescrip(Ev_pos:Ev_pos+3); % F1XX or F2XX (XX is trial such as 01 02 ...)
                
                Tr_pos=[Ev_pos+2:Ev_pos+3];      % Position of Trial
                Trial=str2double(Bdescrip(Tr_pos)); % Trial

                reg_pos=strfind(Bdescrip, 'bf');       % Position of Regresor info
                regType=Bdescrip(reg_pos:reg_pos+1); % Regressor type

                BetaInfo(bi).Block=Block;
                BetaInfo(bi).Trial=Trial;
                BetaInfo(bi).EventType=EventName;
                BetaInfo(bi).regType=regType;
                
            else % from Beta_0002, use the information acquired from Beta_0001
                Bidx=Bdescrip(Bi_pos:Bi_pos+3);                       % 000X
                EventType=Bdescrip(Ev_pos:Ev_pos+1);            % F1 or F2
                if strcmp(EventType,'F1') || strcmp(EventType,'F2')
                    Trial=str2double(Bdescrip(Tr_pos));                % Trial
                    regType=Bdescrip(reg_pos:reg_pos+1);        % Regressor type
                    EventName=Bdescrip(Ev_pos:Ev_pos+3);     % F1XX or F2XX (XX is trial such as 01 02 ...)
                    if Trial==1 && strcmp(EventType,'F1') && strcmp(regType,'bf')
                        Block=Block+1;                                          % Next Block
                    end
                    BetaInfo(bi).Block=Block;
                    BetaInfo(bi).Trial=Trial;
                    BetaInfo(bi).EventType=EventName;
                    BetaInfo(bi).regType=regType;
                end % if EventType
            end %if bi~=1
        end %for bi
    %end %for sub

    clear BetaInfoPM BetaPMNum
    j=1; % to count the number of Beta modulated by [PM_interests]
    for i=1:size(BetaInfo,2)
        if strcmp(BetaInfo(i).regType,PM_interest)       % if the current beta file encodes the activity moulated by the the [PM_interest]
            BetaInfoPM(j)=BetaInfo(i);
            BetaPMNum(j)=i;
            j=j+1;
        else                                                                     % motion or other regressors
        end
    end

    for j0=1:j-1
        BetaInfoPM(j0).BetaPMNum=BetaPMNum(j0); % Add Beta file index (Which Beta files includes the activity modeled by the [PM_interest])
    end
    BetaFilepath=fullfile(betapath, 'BetaInfo.mat');
    save(fullfile(betapath, 'BetaInfo.mat'), 'BetaInfoPM', 'BetaPMNum');
    fprintf ('Beta Labels are saved in %s \n\n', BetaFilepath);
end
