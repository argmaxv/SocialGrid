function [fmri_spec, condcounter] = Model_Grid_PhixB(subj, nblocks, session1, session2, fmri_spec, periodicity)

    [ProjSet, fs, info, ~, ~]=Call_default_PS;
    data_path =ProjSet.DATApath;
 
    for sess = 1:nblocks
        
        conditions = {'F0', 'F1F2', 'Dec', 'Btn'}; % GLM2+GLM3
%         conditions = {'F0', 'F1F2', 'Btn'}; % GLM2
%         conditions = {'F0', 'F1F2', 'Dec', 'Btn'}; % GLM3
        
        cond.conditions = conditions;
        clear onsetlist newses
        if sess<(nblocks/2)+1
            session=session1;
            newses=sess;
        else
            session=session2;
            newses=sess-(nblocks/2);
        end
        condcounter = 0;

%% Specify regressors
        for ev=1:length(conditions)
            condcounter = condcounter+1;    
            fmri_spec.sess(sess).cond(1,condcounter).tmod = 0;

            % Parametric regressors
            Fc=session(newses).F1.correct==1;

            if strcmp(conditions{ev}, 'F0')==1 % F0 onsets and missing trials
                fmri_spec.sess(sess).cond(1,condcounter).onset = [session(newses).F0a.ons; session(newses).F0b.ons; session(newses).F1.ons(~Fc); session(newses).F2.ons(~Fc)];
                fmri_spec.sess(sess).cond(1,condcounter).pmod                            = struct('name', {}, 'param', {}, 'poly', {});
                
            elseif strcmp(conditions{ev}, 'F1F2')==1 % At the time of F1 and F2 onsets
                fmri_spec.sess(sess).cond(1,condcounter).onset = [session(newses).F1.ons(Fc); session(newses).F2.ons(Fc)]; %Correct trials
%% GLM2+GLM3
                idpm=1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)                 = struct('name', 'vF0102', 'param', [cos(periodicity* [ session(newses).F1.pm(Fc,8); session(newses).F2.pm(Fc,7) ])], 'poly', 1); %cos(periodicity*(theta-phi))
                idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'GP', 'param', [session(newses).F1.pm(:,3); session(newses).F2.pm(:,3)], 'poly', 1); % Decision value, GP (growth potential)
                idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'Euc', 'param', [session(newses).F1.pm(Fc,9); session(newses).F2.pm(Fc,8)], 'poly', 1);  %Euclidean distance of the vector
%% GLM2
%                 idpm=1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)                 = struct('name', 'vF0102', 'param', [cos(periodicity* [ session(newses).F1.pm(Fc,8); session(newses).F2.pm(Fc,7) ])], 'poly', 1); %cos(periodicity*(theta-phi))
%% GLM3
%                 idpm=1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)                 = struct('name', 'GP', 'param', [session(newses).F1.pm(:,3); session(newses).F2.pm(:,3)], 'poly', 1); % Decision value, GP (growth potential)
            
            elseif strcmp(conditions{ev}, 'Dec')==1 % At the time of decision (F2 onset)
                fmri_spec.sess(sess).cond(1,condcounter).onset = session(newses).F2.ons(Fc); 
                idpm=1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)                 = struct('name', 'Easiness', 'param', session(newses).F2.pm(Fc,4), 'poly', 1);  %Easiness of the decision |GP1-GP2|

            else % Btn %At the time of button press (motion)
                fmri_spec.sess(sess).cond(1,condcounter).onset = session(newses).(conditions{ev}).ons;
                fmri_spec.sess(sess).cond(1,condcounter).pmod = struct('name', {}, 'param', {}, 'poly', {}); % No parametric modulator, otherwise
            end

            % Duration
            if strcmp(conditions{ev}, 'Btn')==1
                fmri_spec.sess(sess).cond(1,condcounter).duration    = 0; %Stick function
            else 
                fmri_spec.sess(sess).cond(1,condcounter).duration    = 5; %Boxcar function
            end
                fmri_spec.sess(sess).cond(1,condcounter).orth           = 0;  % Default in SPM12 is yes which is 1 (orthoginalized)
        end
        
        for c=1:condcounter
            fmri_spec.sess(sess).cond(c).name = conditions{c};
        end
        
    end
clear session;

%% Add motion regressors
    for sess = 1:nblocks
        fmri_spec.sess(sess).multi           = {''};
        fmri_spec.sess(sess).regress         = struct('name', {}, 'val', {});
        curfd=pwd;
        
        if sess<(nblocks/2)+1
            newses=sess;
            motionregpath=[data_path, info.prefix.day1, subj, fs, info.prefix.Run, num2str(newses)];
            rpfile=dir(fullfile(motionregpath, 'rp*.txt'));
            fmri_spec.sess(sess).multi_reg      = {fullfile(motionregpath, rpfile.name)};  
        else
            newses=sess-(nblocks/2);
            motionregpath=[data_path, info.prefix.day2, subj, fs, info.prefix.Run, num2str(newses)];
            rpfile=dir(fullfile(motionregpath, 'rp*.txt'));
            fmri_spec.sess(sess).multi_reg      = {fullfile(motionregpath, rpfile.name)};          
        end

        fmri_spec.sess(sess).hpf            = 128;
        fmri_spec.fact                      = struct('name', {}, 'levels', {});
        fmri_spec.bases.hrf.derivs          = [0 0];
        fmri_spec.global                    = 'None';
        fmri_spec.cvi                       = 'AR(1)';
        cd(curfd);
    end
end