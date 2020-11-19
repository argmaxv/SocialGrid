function [fmri_spec, condcounter] = Model_NoPhi_5s(subj, nblocks, session1, session2, fmri_spec, periodicity)
    
    [ProjSet, fs, info, ~, ~]=Call_default_PS;
    data_path =ProjSet.DATApath;
    
    for sess = 1:nblocks
        
        conditions = {'F0', 'F1F2', 'Btn'};
        cond.conditions = conditions;
        clear newses
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
            fmri_spec.sess(sess).cond(1,condcounter).tmod              = 0;

            if strcmp(conditions{ev}, 'F0')==1 %F0 and wrong trials
                fmri_spec.sess(sess).cond(1,condcounter).onset             = [session(newses).F0a.ons; session(newses).F0b.ons; session(newses).F1.ons(session(newses).F1.correct==0); session(newses).F2.ons(session(newses).F2.correct==0)];
                fmri_spec.sess(sess).cond(1,condcounter).pmod             = struct('name', {}, 'param', {}, 'poly', {});
                
            elseif strcmp(conditions{ev}, 'F1F2')==1
                % Onset
                fmri_spec.sess(sess).cond(1,condcounter).onset             = [session(newses).F1.ons(session(newses).F1.correct==1); session(newses).F2.ons(session(newses).F2.correct==1)];
                % Parametric regressors
                cosPThF01F02=cos(periodicity*[session(newses).F1.pm(session(newses).F1.correct==1,7); session(newses).F2.pm(session(newses).F2.correct==1,7)]); %Cos periodicity theta
                sinPThF01F02=sin(periodicity*[session(newses).F1.pm(session(newses).F1.correct==1,7); session(newses).F2.pm(session(newses).F2.correct==1,7)]);%Sin periodicity theta
                idpm=1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)                    = struct('name', 'CvF01F02', 'param', cosPThF01F02, 'poly', 1); %Cos vector
                idpm=idpm+1; fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)         = struct('name', 'SvF01F02', 'param', sinPThF01F02 , 'poly', 1); %Sin vector
             
            else % Btn - Stick function with no PM
                fmri_spec.sess(sess).cond(1,condcounter).onset             = session(newses).(conditions{ev}).ons;
                fmri_spec.sess(sess).cond(1,condcounter).pmod             = struct('name', {}, 'param', {}, 'poly', {});
            end

            % Duration            
            if strcmp(conditions{ev}, 'F1F2')==1
                fmri_spec.sess(sess).cond(1,condcounter).duration    = 5; % Boxcar, Presentation duration of F1 and F2
            else % F0 and Btn
                fmri_spec.sess(sess).cond(1,condcounter).duration    = 0; % Stick function
            end
            
            fmri_spec.sess(sess).cond(1,condcounter).orth	= 0; % Default in SPM12 is yes which is 1 (orthoginalized)
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