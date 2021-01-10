%% Grid Category
    % The inferred vectors (Theta-Phi) are categorized to the bins
    % according to the periodicity of interests
    % It returns AltPeriodicity which includes 1) periodicities :4~8
    % 2) onoff : 1, aligned (on-grid); 0, misaligned (off-grid)
    % 3) cat: category of the bin (total num. bin = periodicity*2-1 so, 11 bins in six-fold periodicity)

clear; close all; clc
[ProjSet, fs, info, ROI, fname]=Call_default_PS;
%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%
phiset.ROI = ROI.Grid{1};                       % which grid angle will extract (e.g. EC_Grid)
phiset.type = 'PhixB';                              % cross validation of Phi which is estimated differnt blocks (across blocks);
%phiset.type = 'PhixB';                           % Alternatively, PhixD - across days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AltPeriodicity]=Altanative_Periodicity(phiset);
svpath=[ProjSet.PhiInfopath, phiset.ROI, fs];
gridcatpath=fullfile(svPath,  ['GridCat_', phiset.ROI, '_', phiset.type,'.mat']);
save(gridcatpath, 'AltPeriodicity'); %  e.g. GridCat_EC_Grid_PhixB.mat

function  [AltPeriodicity]=Altanative_Periodicity(phiset)
    phiset1=phiset.ROI;
    phiset2=phiset.type;
    [ProjSet, fs, info, ~, ~]=Call_default_PS;
    [sub0, subn] = CallSubj_PS;
    onsetpath=[ProjSet.Metapath, phiset1, fs, phiset2, fs]; % Where Onsets are (behavioral data)
    periodicity=4:8; % including the alternative periodicity
    nBlock=info.Nses*info.Nday; % total num block
    svPath=[ProjSet.PhiInfopath, phiset1, fs];
    
    for sbj=1:subn
        Day1mat=dir(fullfile(onsetpath, ['*', info.prefix.day1, sub0{sbj}, '.mat']));
        Day2mat=dir(fullfile(onsetpath, ['*', info.prefix.day2, sub0{sbj}, '.mat']));

        RadVec=[];
        for ses=1:nBlock
            clear cur_ses curOnset;
            if ses<=info.Nses % if data was acquired from Day1 scan
                cur_ses=ses;
                curOnset=load(fullfile(Day1mat.folder, Day1mat.name));
            else                       % if data was acquired from Day2 scan
                cur_ses=ses-info.Nses;
                curOnset=load(fullfile(Day2mat.folder, Day2mat.name));
            end

            vec1=find(contains(curOnset.Session.se(cur_ses).F1.pmidx, 'RPhiF01'));  
            vec2=find(contains(curOnset.Session.se(cur_ses).F2.pmidx, 'RPhiF02'));  
            RadVec1=curOnset.Session.se(cur_ses).F1.pm(:, vec1);    %Rad angle Theta-Phi F0F1
            RadVec2=curOnset.Session.se(cur_ses).F2.pm(:, vec2);    %Rad angle Theta-Phi F0F2
            RadVec0=[RadVec1; RadVec2];                   % Match the format of the Beta files (each block)
            RadVec=[RadVec; RadVec0];                       % Match the format of the Beta files (each subject)        

        end
            TheRad(:,sbj)=RadVec;                                 % Match the format of the Beta files (NumVecRad x subjects)
            for pr=1:numel(periodicity) %=4:8
                [on, cat]=isRGridN(RadVec,periodicity(pr));
                AltPeriodicity(pr).periodicity=periodicity(pr);
                AltPeriodicity(pr).onoff(:,sbj)=on;
                AltPeriodicity(pr).cat(:,sbj)=cat;
            end
    end
end
