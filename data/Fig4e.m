close all;
clear all;
clc
load('GP_Beta_extracted.mat');
GP_Y=[GPon_Y,GPoff_Y];
figure; bar([mean(GPon_Y), mean(GPoff_Y)]);
hold on;
errorbar([1:6],[mean(GPon_Y), mean(GPoff_Y)], [SEM(GPon_Y), SEM(GPoff_Y)], 'o');
scaatX=repmat([1:6],21,1);
hold on;
scatter(scaatX(:),GP_Y(:),'.');
xlim([0 7]); ylim([.0 .6]);
xticklabels({'dmPFC','rTPJ','lTPJ','dmPFC','rTPJ','lTPJ'});

figure; GPdiff=GPon_Y-GPoff_Y; bar([mean(GPdiff)]);
hold on;
errorbar([1:3], mean(GPdiff), SEM(GPdiff), '.');
[h,p,ci,stat]=ttest(GPdiff);
[h,pfwe]=bonferroni_holm(p);
pfwe
stat.tstat
