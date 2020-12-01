
function [val, meanval] = extract_meanactivity_fromROI(mask,data)
% from *_roi.mat
    V = spm_vol(deblank(data));
    roi = maroi(deblank(mask));
    xyz = voxpts(roi,deblank(data));
    val = spm_sample_vol(V,xyz(1,:),xyz(2,:),xyz(3,:),0);
    meanval=nanmean(val);
end