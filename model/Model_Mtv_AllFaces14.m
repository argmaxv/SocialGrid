function [fmri_spec, condcounter] = Model_Mtv_AllFaces14(subj, nblocks, session, fmri_spec)

    [ProjSet, fs, info, ~, ~]=Call_default_PS;
    data_path = ProjSet.DATApath;
    
    for sess = 1:nblocks
%         conditions = {'Face1', 'Face2', 'Face3','Face4','Face5','Face6','Face7','Face8','Face9','Face10', 'Face11','Face12','Face13','Face14','Face15','Face16'};
        conditions = {'Face2', 'Face3','Face4','Face5','Face6','Face7','Face8','Face9','Face10', 'Face11','Face12','Face13','Face14','Face15'};
        cond.conditions = conditions;
        onsetlist=fieldnames(session(sess));
        condcounter = 0;   
        
%% Specify regressors
        for ev=1:length(conditions)         
            condcounter = condcounter+1;  

% Onset
            fmri_spec.sess(sess).cond(1,condcounter).onset             = session(sess).Fons.(conditions{ev}).ons;
            fmri_spec.sess(sess).cond(1,condcounter).tmod              = 0;
            fmri_spec.sess(sess).cond(1,condcounter).pmod           = struct('name', {}, 'param', {}, 'poly', {});  
   
% Parametric regressors
        for n=2:15
            if strcmp(conditions{ev}, ['Face' num2str(n)])==1
                curFace=['Face' num2str(n)];
                fmri_spec.sess(sess).cond(1,condcounter).pmod                                         = struct('name', {}, 'param', {}, 'poly', {}); %
            end
        end  
        
% Duration    
                fmri_spec.sess(sess).cond(1,condcounter).duration    = 2.5; %F0 was presented for 2.5 sec and F1 and F2 were presented for 5 sec. so, duration was set to 2.5s.

                fmri_spec.sess(sess).cond(1,condcounter).orth                = 0; %default in SPM12 is yes which is 1 (orthoginalized)
                for c=1:condcounter
                fmri_spec.sess(sess).cond(1,c).name = conditions{c};
                end
        end %ev
    end % for session blocks 
        clear session;
%% Add motion regressors
        for sess = 1:nblocks
            fmri_spec.sess(sess).multi           = {''};
            fmri_spec.sess(sess).regress         = struct('name', {}, 'val', {});
            motionregpath=[data_path, subj, fs, info.prefix.Run, num2str(sess)];
            rpfile=dir(fullfile(motionregpath, 'rp*.txt'));
            fmri_spec.sess(sess).multi_reg      = {fullfile(motionregpath, rpfile.name)}; 
            fmri_spec.sess(sess).hpf            = 128;
            fmri_spec.fact                      = struct('name', {}, 'levels', {});
            fmri_spec.bases.hrf.derivs          = [0 0];
            fmri_spec.global                    = 'None';
            fmri_spec.cvi                       = 'AR(1)';
            cd(curfd);
        end
    end