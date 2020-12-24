%% RSA Searchlight 
% 0. Set n serchlight including 100 voxels
% by SPARK 10.Oct.2018

[ProjSet, fs, info, ROI, fname]=Call_default_PS;
maskpath=ProjSet.Maskpath;
Vmask=spm_vol(fullfile(maskpath, fname.maskname));
Vmask.data=spm_read_vols(Vmask);
Ll = rsa.defineSearchlight({Vmask},Vmask,'shere',[15 100]);
save(fullfile(maskpath,'searchlight_100.mat') ,'-struct', 'Ll' );