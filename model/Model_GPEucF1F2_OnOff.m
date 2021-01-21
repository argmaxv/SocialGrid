function [fmri_spec, condcounter] = Model_GPEucF1F2_OnOff(subj, nblocks, session1, session2, fmri_spec, OnOff_sub)

    [ProjSet, fs, info, ~, ~]=Call_default_PS;
    data_path = ProjSet.DATApath;
     
    for sess = 1:nblocks
        clear On_Off
        On_Off=OnOff_sub.block(sess).OnOff;
        conditions = {'F12on', 'F12off', 'Btn'};
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
            if strcmp(conditions{ev}, 'F12on')==1
                clear temponset tempGP
                temponset=[session(newses).F1.ons; session(newses).F2.ons];
                tempGP=[session(newses).F1.pm(:,3); session(newses).F2.pm(:,3)];
                tempEuc=[session(newses).F1.pm(:,9); session(newses).F2.pm(:,8)];
                fmri_spec.sess(sess).cond(1,condcounter).onset = temponset(On_Off);
                fmri_spec.sess(sess).cond(1,condcounter).pmod                            = struct('name', {}, 'param', {}, 'poly', {});
                idpm=1;             fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'GPon', 'param', tempGP(On_Off), 'poly', 1);
                idpm=idpm+1;   fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'Eucon', 'param', tempEuc(On_Off), 'poly', 1);
                fmri_spec.sess(sess).cond(1,condcounter).duration    = 5;
            elseif strcmp(conditions{ev}, 'F12off')==1
                clear temponset tempGP
                temponset=[session(newses).F1.ons; session(newses).F2.ons];
                tempGP=[session(newses).F1.pm(:,3); session(newses).F2.pm(:,3)];
                tempEuc=[session(newses).F1.pm(:,9); session(newses).F2.pm(:,8)];
                fmri_spec.sess(sess).cond(1,condcounter).onset = temponset(~On_Off);
                fmri_spec.sess(sess).cond(1,condcounter).pmod                            = struct('name', {}, 'param', {}, 'poly', {});
                idpm=1;             fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'GPoff', 'param', tempGP(~On_Off), 'poly', 1);
                idpm=idpm+1;   fmri_spec.sess(sess).cond(1,condcounter).pmod(idpm)      = struct('name', 'Eucoff', 'param', tempEuc(~On_Off), 'poly', 1);
                fmri_spec.sess(sess).cond(1,condcounter).duration    = 5;
           
            else % F0, Btn
                fmri_spec.sess(sess).cond(1,condcounter).onset = session(newses).(conditions{ev}).ons;
                fmri_spec.sess(sess).cond(1,condcounter).pmod = struct('name', {}, 'param', {}, 'poly', {}); % No parametric modulator, otherwise
                fmri_spec.sess(sess).cond(1,condcounter).duration    = 0;
            end
            
            fmri_spec.sess(sess).cond(1,condcounter).orth       = 0; %default in SPM12 is yes which is 1 (orthoginalized)
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