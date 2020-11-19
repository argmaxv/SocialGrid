
function zscore_nii(path, file, prefix)

% Getting the input of the 1st level contrast, it returns the zscored map in
% the same folder.

        curnii=fullfile(path, file);
        mask=spm_read_vols(spm_vol( fullfile(path, 'mask.nii')));
        struct=spm_vol(curnii);
        vols=spm_read_vols(struct);
        vols(mask~=1)=NaN;        
        zvols=(vols-mean(vols(:),'omitnan'))/std(vols(:),'omitnan');
        struct.fname=fullfile(path, [prefix, file(4:end)]);
        spm_write_vol(struct, zvols);
end