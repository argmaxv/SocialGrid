function x=MultiErrorbar(model_series)
% -----------------------------------------------------
% Compute the x axis values to ground errorbars on multiple bars
%
% e.g.
% bar(RtMean); hold on
% x=MultiErrorbar(RtMean);
% errorbar(x,RtMean, RtSE, 'k', 'linestyle', 'none');
%
% -----------------------------------------------------

    ngroups = size(model_series, 1);
    nbars = size(model_series, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        for j=1:ngroups
            x(j,i)=j- groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        end
    end
end