%% RSA searchlight 
% 4. 2nd level analysis

close all
clear
clc

[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[Subject, NSubject]=CallSubj_PS; 
ModelNames={'Mtv_FaceNMatch14', 'Mtv_AllFaces14'};

    for Mi=1:numel(ModelNames)
        sL_model=ModelNames{Mi};
        scanpath=[ProjSet.Respath, sL_model, fs, 'sL', fs, 'TauMap_Euc_T'];
        lv2path{1}=[scanpath, fs, 'Lv2/'];
        lv2conpath=[scanpath, fs, 'Lv2', fs, 'con'];
        if ~exist(lv2path{1})
            mkdir(lv2path{1});
            mkdir(lv2conpath);
        end
        factorial_design.dir = lv2path;
        smoothedmap=[info.prefix.smooth, sL_model, '_*.nii'];
        list=dir(fullfile(scanpath, smoothedmap));
        
       sia=0;
        for si=1:numel(list)
                sia=sia+1;
                factorial_design.des.t1.scans{sia,1}=strcat(list(si).folder, filesep, list(si).name);
                copyfile( fullfile(list(si).folder, list(si).name), strcat(lv2conpath, fs, sprintf('con_%02d.nii',sia)) ,'f');
        end
        factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        factorial_design.masking.tm.tm_none = 1;
        factorial_design.masking.im = 1;
        factorial_design.masking.em = {''};
        factorial_design.globalc.g_omit = 1;
        factorial_design.globalm.gmsca.gmsca_no = 1;
        factorial_design.globalm.glonorm = 1;
        matlabbatch{1}.spm.stats.factorial_design=factorial_design;
        
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        
        
        con.spmmat = {[factorial_design.dir{1}, fs, 'SPM.mat']};
        tconname={['RDM_', Models{Mi}]};
        con.consess{1}.tcon.name = tconname{1};
        con.consess{1}.tcon.weights = 1;
        con.consess{1}.tcon.sessrep = 'none';
        con.delete = 0;
        matlabbatch{3}.spm.stats.con=con;

        results.spmmat = con.spmmat ;
        for tconN=1:numel(tconname)
            results.conspec(tconN).titlestr = '';
            results.conspec(tconN).contrasts = tconN;
            results.conspec(tconN).threshdesc = 'none';
            results.conspec(tconN).thresh = 0.001;
            results.conspec(tconN).extent = 0;
            results.conspec(tconN).conjunction = 1;
            results.conspec(tconN).mask.none = 1;
        end
        results.units = 1;
        results.export{1}.ps = false;
        matlabbatch{4}.spm.stats.results=results;

        spm('defaults', 'FMRI');
        spm_jobman('run',matlabbatch);
        clear matlabbatch;
    end
