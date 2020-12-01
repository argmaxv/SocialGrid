%% Phi_histogram
% Perfomming Rayleigh tests and drawing distributon of Phi estimated across voxels

clear; close all; clc

[ProjSet, fs, info, ROI, fname]=Call_default_PS;
[subj, subn] = CallSubj_PS;
Phipath=ProjSet.PhiInfopath;
ROIs=ROI.Grid;
periodicity=info.periodicity;
design_name = ['Grid', num2str(periodicity), '_F01F02_5s']; %Make the same with it defined in Get_Phi.m
Nday=info.Nday;
Nses=info.Nses;

for roix=1:1%numel(ROIs)
    clear betaname phi_roi meanPhi_indv0 meanPhi_indv pval z pval_day
    
    % load results of Get-Phi.m
    Phifilename = ['Phi_', design_name, '_', ROIs{roix}, '.mat'];
    betaname=dir(fullfile([Phipath, fs, ROIs{roix}], Phifilename));
    load(fullfile(betaname.folder, betaname.name));
    PhiInfo(roix).ROI=ROIs{roix};
    eval(['Lv1Model=', design_name, ';'])
    
    for s=1:subn
	% Compute the circular mean across sessions per subject
    % Note that raw 'Phi' values are not yet devided by periodicity.
    % That is, Rayleigh test is performed with these raw Phi but the grid
    % angle should be computed by deviding six-fold periodicity
        clear day
        for dayi=1:Nday
            day(dayi).Phi=[];
        end
        for bl=1:Nses*Nday
            if bl<=Nses
                day(1).Phi=[day(1).Phi; Lv1Model(s).vF01F02.ses(bl).Phi];
            else
                day(2).Phi=[day(2).Phi; Lv1Model(s).vF01F02.ses(bl).Phi];
            end
        end

        for dayi=1:Nday
            PhiInfo(roix).meanPhi(s,dayi)=circ_mean(circ_mean(day(dayi).Phi)');
            [pval_day(s,dayi), zval_day(s,dayi)] = circ_rtest(circ_mean(day(dayi).Phi)'); % Rayleigh test across voxels while Phi of each voxel is the sirular mean of Phi across blocks
            [pval_voxses(s,dayi), zval_voxses(s,dayi)] = circ_rtest(day(dayi).Phi(:));       % Rayleigh test across voxels in all blocks
            meanPhi{s,dayi}=circ_mean(day(dayi).Phi)';
            allPhi{s,dayi}=day(dayi).Phi(:);
        end
        
    end %for subj s
    
    % PolarHistogram
    % Note that the Phi range (angle) should be  [0, pi/12, 2*pi/12,
    % 3*pi/12] not [0, pi/2, 2*pi/2, 3*pi/2], since Phi was not devided by
    % periodicity on here.
    for dayi=1:Nday
        Phifig((roix-1)*Nday+dayi)=figure; 
        for s=1:subn
            subplot(5,5,s);
            polarhistogram(allPhi{s,dayi});                % voxels x blocks
            %polarhistogram(meanPhi{s,dayi});       % voxels (circular mean
            %across blocks)
            phidist=gca;
            thetaticks(0:90:270);
            %phidist.ThetaTickLabel=[0, pi/12, 2*pi/12, 3*pi/12];
            phidist.ThetaAxisUnits='radians';
        end
        suptitle([PhiInfo(roix).ROI(1:end-4), ' Day ', num2str(dayi)]);
    end
        
end