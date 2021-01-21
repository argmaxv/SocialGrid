function [fmri_spec, condcounter] = Model_EachDec_GP(subj, nblocks, session1, session2, fmri_spec)
 
    [ProjSet, fs, info, ~, ~]=Call_default_PS;
    data_path = ProjSet.DATApath;
    
    for sess = 1:nblocks
        
        conditions = {'F1', 'F2', 'Btn'};
        cond.conditions = conditions;
        clear newses
        if sess<(nblocks/2)+1
            session=session1;
            newses=sess;
        else
            session=session2;
            newses=sess-(nblocks/2);
        end
        condcounter1= 0;
        condcounter2= 0;

%% Specify regressors
        for ev=1:length(conditions)
            if strcmp(conditions{ev}, 'F1')==1
                for trial=1:length(session(newses).(conditions{ev}).ons)
                	condcounter1 = condcounter1+1;
                    fmri_spec.sess(sess).cond(condcounter1).name          = [conditions{ev}, sprintf('%02d', condcounter1)];
                    fmri_spec.sess(sess).cond(1,condcounter1).onset       = session(newses).(conditions{ev}).ons(trial); 
                    fmri_spec.sess(sess).cond(1,condcounter1).tmod        = 0;
                    fmri_spec.sess(sess).cond(1,condcounter1).pmod       = struct('name', {}, 'param', {}, 'poly', {});
                    idpm=1;
                    fmri_spec.sess(sess).cond(1,condcounter1).pmod(idpm)       = struct('name', 'GPF101', 'param', session(newses).(conditions{ev}).pm(trial,3) , 'poly', 1);
                    fmri_spec.sess(sess).cond(1,condcounter1).duration    = 5;
                    fmri_spec.sess(sess).cond(1,condcounter1).orth           = 0; %default in SPM12 is yes which is 1 (orthoginalized)
                end
            elseif strcmp(conditions{ev}, 'F2')==1
                
                for trial=1:length(session(newses).(conditions{ev}).ons)
                	condcounter2 =condcounter2+1;
                    fmri_spec.sess(sess).cond(condcounter1+condcounter2).name           = [conditions{ev}, sprintf('%02d', condcounter2)];
                    fmri_spec.sess(sess).cond(1,condcounter1+condcounter2).onset        = session(newses).(conditions{ev}).ons(trial); 
                    fmri_spec.sess(sess).cond(1,condcounter1+condcounter2).tmod         = 0;
                    fmri_spec.sess(sess).cond(1,condcounter1+condcounter2).pmod        = struct('name', {}, 'param', {}, 'poly', {});
                    idpm=1;
                    fmri_spec.sess(sess).cond(1,condcounter1+condcounter2).pmod(idpm)    = struct('name', 'GPF202', 'param', session(newses).(conditions{ev}).pm(trial,3) , 'poly', 1);
                    fmri_spec.sess(sess).cond(1,condcounter1+condcounter2).duration    = 5;
                    fmri_spec.sess(sess).cond(1,condcounter1+condcounter2).orth                = 0; %default in SPM12 is yes which is 1 (orthoginalized)
                end                    
            else %Btn
                condcounter=1+condcounter1+condcounter2;
                fmri_spec.sess(sess).cond(condcounter).name                = conditions{ev};
                fmri_spec.sess(sess).cond(1,condcounter).onset             = session(newses).(conditions{ev}).ons; 
                fmri_spec.sess(sess).cond(1,condcounter).tmod              = 0;
                fmri_spec.sess(sess).cond(1,condcounter).pmod             = struct('name', {}, 'param', {}, 'poly', {});
                fmri_spec.sess(sess).cond(1,condcounter).duration         = 0;
                fmri_spec.sess(sess).cond(1,condcounter).orth                = 0; %default in SPM12 is yes which is 1 (orthoginalized)
            end
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
            cd([data_path,  info.prefix.day1, subj, fs, info.prefix.Run, num2str(newses)]);
            rpfile=dir('rp*.txt');
            fmri_spec.sess(sess).multi_reg	= {[data_path, info.prefix.day1, subj, fs, info.prefix.Run, num2str(newses), fs rpfile.name]}; 
        else
            newses=sess-(nblocks/2);
            cd([data_path,  info.prefix.day2, subj, fs, info.prefix.Run, num2str(newses)]);
            rpfile=dir('rp*.txt');
            fmri_spec.sess(sess).multi_reg	= {[data_path, info.prefix.day2, subj, fs, info.prefix.Run, num2str(newses), fs rpfile.name]};
        end

        fmri_spec.sess(sess).hpf            = 128;
        fmri_spec.fact                      = struct('name', {}, 'levels', {});
        fmri_spec.bases.hrf.derivs          = [0 0];
        fmri_spec.global                    = 'None';
        fmri_spec.cvi                       = 'AR(1)';
        cd(curfd);
    end
    
end